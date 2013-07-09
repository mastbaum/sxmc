#include <iostream>
#include "pdfz.h"
#include <math.h>
#include <cuda.h>
#include <math_constants.h> // CUDA header
#include <TStopwatch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>

#include "cuda_compat.h"

namespace pdfz {
    const int MAX_NFIELDS = 10;

    // Changes the size of an array while preserving its contents.
    // Shrinking the array only preserves the initial entries, while truncating
    // those past the new length.
    template <typename T>
    void resize_array(hemi::Array<T> &array, const size_t n)
    {
        std::vector<T> tmp(n);
        T *orig_contents = array.hostPtr();

        for (unsigned int i=0; i < array.size() && i < n; i++)
            tmp[i] = orig_contents[i];

        array.copyFromHost(&tmp.front(), n); // Resizes array automatically
    }

    // Hidden CUDA state for evaluator.  Not visible to Eval class users.
    struct CudaState
    {
        cudaStream_t stream;
    };


    struct SystematicDescriptor
    {
        short type;
        short obs, extra_field;
        short par;
    };

    Eval::Eval(const std::vector<float> &_samples, int _nfields, int _nobservables,
               const std::vector<float> &_lower, const std::vector<float> &_upper) :
        nfields(_nfields), nobservables(_nobservables),
        lower(_lower.size(), true), upper(_upper.size(), true), syst(0)
    {
        if (_samples.size() % _nfields != 0)
            throw Error("Length of samples array is not divisible by number of fields.");

        if (_nobservables == 0 )
            throw Error("Number of observables in PDF is zero.");

        if (_nobservables > nfields )
            throw Error("Number of observables cannot be greater than number of fields.");

        if ( (int) _upper.size() != _nobservables)
            throw Error("Number of upper bounds must be same as number of observables.");

        if ( (int) _lower.size() != _nobservables)
            throw Error("Number of lower bounds must be same as number of observables.");

        // Transfer lower and upper bounds to hemi arrays for use on device
        this->lower.copyFromHost(&_lower.front(), _lower.size());
        this->upper.copyFromHost(&_upper.front(), _upper.size());

        // Setup the CUDA state class
        this->cuda_state = new CudaState;
        #ifdef __CUDACC__
        checkCuda( cudaStreamCreate(&(this->cuda_state->stream)) );
        #else
        this->cuda_state->stream = 0;
        #endif

        // The handling of samples is up to the Eval subclass!
    }

    Eval::~Eval()
    {
        delete this->cuda_state;
        delete this->syst;
    }

    void Eval::SetPDFValueBuffer(hemi::Array<float> *output, int offset, int stride)
    {
        this->pdf_buffer = output;
        this->pdf_offset = offset;
        this->pdf_stride = stride;
    }

    void Eval::SetNormalizationBuffer(hemi::Array<unsigned int> *norm, int offset)
    {
        this->norm_buffer = norm;
        this->norm_offset = offset;
    }

    void Eval::SetParameterBuffer(hemi::Array<float> *params, int offset, int stride)
    {
        this->param_buffer = params;
        this->param_offset = offset;
        this->param_stride = stride;
    }

    void Eval::AddSystematic(const Systematic &syst)
    {
        if (this->syst)
            resize_array(*(this->syst), this->syst->size() + 1);
        else
            this->syst = new hemi::Array<SystematicDescriptor>(1, true);

        SystematicDescriptor desc;
        desc.type = syst.type;
        if (syst.type == Systematic::SHIFT) {
            const ShiftSystematic &shift = dynamic_cast<const ShiftSystematic &>(syst);
            desc.obs = shift.obs;
            desc.par = shift.par;
        } else if (syst.type == Systematic::SCALE) {
            const ScaleSystematic &scale = dynamic_cast<const ScaleSystematic &>(syst);
            desc.obs = scale.obs;
            desc.par = scale.par;
        } else if (syst.type == Systematic::RESOLUTION_SCALE) {
            const ResolutionScaleSystematic &res = dynamic_cast<const ResolutionScaleSystematic &>(syst);
            desc.obs = res.obs;
            desc.extra_field = res.true_obs;
            desc.par = res.par;
        } else {
            throw Error("Unknown systematic type");
        }

        this->syst->writeOnlyHostPtr()[this->syst->size() - 1] = desc;
    }


    ///////////////////// EvalHist ///////////////////////

    EvalHist::EvalHist(const std::vector<float> &_samples, int nfields, int nobservables,
                       const std::vector<float> &lower, const std::vector<float> &upper,
                       const std::vector<int> &_nbins, bool optimize) :
        Eval(_samples, nfields, nobservables, lower, upper),
        samples(_samples.size(), false), read_bins(0), 
        nbins(_nbins.size(), true), bin_stride(_nbins.size(), true), bins(0), needs_optimization(optimize)
    {
        if ( (int) _nbins.size() != nobservables)
            throw Error("Size of nbins array must be same as number of observables.");

        if (nfields > MAX_NFIELDS)
            throw Error("Exceeded maximum number of fields per sample.  Edit MAX_NFIELDS in pdfz.cpp to fix this!");

        this->samples.copyFromHost(&_samples.front(), _samples.size());
        this->nbins.copyFromHost(&_nbins.front(), _nbins.size());

        // Compute bin volume
        this->bin_volume = 1.0f;
        for (int i=0; i < nobservables; i++)
            this->bin_volume *= (upper[i] - lower[i]) / _nbins[i];


        // Compute stride for row-major order storage
        const int ndims = (int) this->bin_stride.size();
        int *bin_stride_host = this->bin_stride.writeOnlyHostPtr();

        bin_stride_host[ndims - 1] = 1;
        for (int i=ndims-2; i >= 0; i--)
            bin_stride_host[i] = _nbins[i+1] * bin_stride_host[i+1];

        this->total_nbins = bin_stride_host[0] * _nbins[0];

        if (this->total_nbins == 0)
            throw Error("Cannot make histogram with zero bins.");

        this->bins = new hemi::Array<unsigned int>(this->total_nbins, true);


        this->bin_nthreads_per_block = 256;
        if (_samples.size() > 100000)
            this->bin_nblocks = 64;
        else if (_samples.size() > 10000)
            this->bin_nblocks = 16;
        else
            this->bin_nblocks = 4;

        this->eval_nthreads_per_block = 256;
        this->eval_nblocks = 64;
    }

    EvalHist::~EvalHist()
    {
        delete this->read_bins;
        delete this->bins;
    }

    void EvalHist::SetEvalPoints(const std::vector<float> &points)
    {
        if (points.size() % this->nobservables != 0)
            throw Error("Number of entries in evaluation points array not divisible by number of observables.");

        delete this->read_bins;
        this->read_bins = new hemi::Array<int>(points.size(), false);
        int *read_bins = this->read_bins->writeOnlyHostPtr();

        // Precompute the bin number corresponding to each evaluation point
        // *** This never changes between PDF evaluations! ***


        // Pointer aliases to Hemi array contents (for convenience)
        const float *lower = this->lower.readOnlyHostPtr();
        const float *upper = this->upper.readOnlyHostPtr();
        const int *nbins = this->nbins.readOnlyHostPtr();
        const int *bin_stride = this->bin_stride.readOnlyHostPtr();

        std::vector<float> bin_scale(this->nobservables);

        for (int iobs=0; iobs < this->nobservables; iobs++) {
            float span = upper[iobs] - lower[iobs];
            bin_scale[iobs] = nbins[iobs] / span;
        }

        for (unsigned int ipoint=0; ipoint < points.size(); ipoint++) {
            bool in_pdf_domain = true;
            int bin_id = 0;

            for (int iobs=0; iobs < this->nobservables; iobs++) {
                int ielement = this->nobservables * ipoint + iobs;
                float element = points[ielement];

                // Throw out this event if outside of PDF domain
                if (element < lower[iobs] || element >= upper[iobs]) {
                    in_pdf_domain = false;
                    break;
                }

                bin_id += (int)( (element - lower[iobs]) * bin_scale[iobs] ) * bin_stride[iobs];
            }

            if (in_pdf_domain)
                read_bins[ipoint] = (int) bin_id;
            else
                read_bins[ipoint] = -1; // Filled in with NaN during evaluation
        }
    }

    ///// EvalHist kernels
    HEMI_DEV_CALLABLE_INLINE
    void apply_systematic(const SystematicDescriptor *syst, float *fields, const float *parameters, const int param_stride)
    {
        switch (syst->type) {
            case Systematic::SHIFT:
            fields[syst->obs] += parameters[syst->par * param_stride];
            break;
            case Systematic::SCALE:
            fields[syst->obs] *= (1 + parameters[syst->par * param_stride]);
            break;
            case Systematic::RESOLUTION_SCALE:
            fields[syst->obs] += parameters[syst->par * param_stride] * (fields[syst->obs] - fields[syst->extra_field]);
            break;
        }
    }

    HEMI_KERNEL(zero_hist)(int total_nbins, unsigned int *bins, unsigned int *norm)
    {
        int offset = hemiGetElementOffset();
        int stride = hemiGetElementStride();

        if (offset == 0)
            *norm = 0;

        for (int i=offset; i < total_nbins; i += stride)
            bins[i] = 0;
    }

    HEMI_KERNEL(bin_samples)(int ndata, const float *data,
                             const int nobs, const int nfields, 
                             const int * __restrict__ bin_stride, const int * __restrict__ nbins,
                             const float * __restrict__ lower, const float * __restrict__ upper,
                             const int nsyst, const SystematicDescriptor * __restrict__ syst,
                             const float * __restrict__ parameters, const int param_stride,
                             unsigned int *bins, unsigned int *norm)
    {
        int offset = hemiGetElementOffset();
        int stride = hemiGetElementStride();
        float field_buffer[MAX_NFIELDS];
        const int nsamples = ndata / nfields;

        float bin_scale[MAX_NFIELDS];
        for (int iobs=0; iobs < nobs; iobs++)
          bin_scale[iobs] = nbins[iobs] / (upper[iobs] - lower[iobs]);

        unsigned int thread_norm = 0;

        for (int isample=offset; isample < nsamples; isample += stride) {
            bool in_pdf_domain = true;
            int bin_id = 0;

            // Copy fields
            for (int ifield=0; ifield < nfields; ifield++)
                field_buffer[ifield] = data[isample * nfields + ifield];

            // Apply systematics
            for (int isyst=0; isyst < nsyst; isyst++)
                apply_systematic(syst + isyst, field_buffer, parameters, param_stride);

            // Compute histogram bin
            for (int iobs=0; iobs < nobs; iobs++) {
                float element = field_buffer[iobs];
                // Throw out this event if outside of PDF domain
                if (element < lower[iobs] || element >= upper[iobs]) {
                    in_pdf_domain = false;
                    break;
                }

                bin_id += (int)( (element - lower[iobs]) * bin_scale[iobs] ) * bin_stride[iobs];
            }

            // Add to histogram if sample in PDF domain
            if (in_pdf_domain) {
                atomicAdd(bins + bin_id, 1);
                thread_norm += 1;
            }
        }

        atomicAdd(norm, thread_norm);
    }

    HEMI_KERNEL(eval_pdf)(int npoints, const int *read_bins,
                          const unsigned int * __restrict__ bins,
                          const unsigned int * __restrict__ norm,
                          float bin_volume,
                          float *output, int output_stride)
    {
        int offset = hemiGetElementOffset();
        int stride = hemiGetElementStride();
        const float bin_norm = *norm * bin_volume;

        for (int ipoint=offset; ipoint < npoints; ipoint += stride) {
            int bin_id = read_bins[ipoint];

            float pdf_value = 0.0f;
            if (bin_id < 0)
                pdf_value = nanf("");
            else
                pdf_value = bins[bin_id] / bin_norm;

            output[output_stride * ipoint] = pdf_value;
        }
    }
    ///// End EvalHist kernels

    void EvalHist::EvalAsync(bool do_eval_pdf)
    {
        if (this->needs_optimization)
            this->Optimize();

        // Handle no systematics case
        int nsyst = 0;
        const SystematicDescriptor *syst_ptr = 0;
        if (this->syst) {
            nsyst = this->syst->size();
            syst_ptr = this->syst->readOnlyPtr();
        }

        HEMI_KERNEL_LAUNCH(zero_hist, this->eval_nblocks, this->eval_nthreads_per_block, 0, this->cuda_state->stream,
                           this->total_nbins, this->bins->writeOnlyPtr(), this->norm_buffer->writeOnlyPtr() + this->norm_offset);
        HEMI_KERNEL_LAUNCH(bin_samples, this->bin_nblocks, this->bin_nthreads_per_block, 0, this->cuda_state->stream,
                           (int) this->samples.size(), this->samples.readOnlyPtr(), 
                           this->nobservables, this->nfields,
                           this->bin_stride.readOnlyPtr(), this->nbins.readOnlyPtr(),
                           this->lower.readOnlyPtr(), this->upper.readOnlyPtr(),
                           nsyst, syst_ptr,
                           this->param_buffer->readOnlyPtr() + this->param_offset, this->param_stride,
                           this->bins->ptr(), this->norm_buffer->writeOnlyPtr() + this->norm_offset);

        if (this->read_bins == 0 || !do_eval_pdf)
            return; // This can happen if someone wants to create a histogram with no eval points.

        HEMI_KERNEL_LAUNCH(eval_pdf, this->eval_nblocks, this->eval_nthreads_per_block, 0, this->cuda_state->stream,
                           (int) this->read_bins->size(), this->read_bins->readOnlyPtr(),
                           this->bins->readOnlyPtr(), this->norm_buffer->readOnlyPtr() + this->norm_offset,
                           this->bin_volume,
                           this->pdf_buffer->writeOnlyPtr() + this->pdf_offset, this->pdf_stride);
    }

    void EvalHist::EvalFinished()
    {
        #ifdef __CUDACC__
        checkCuda( cudaStreamSynchronize(this->cuda_state->stream) );
        #endif
    }

    TH1* EvalHist::CreateHistogram()
    {
        if (this->nobservables > 3)
            throw Error("Cannot EvalHist::CreateHistogram for dimensions greater than 3!");

        bool orig_optimization_flag = this->needs_optimization;
        this->needs_optimization = false; // never optimize when making histogram!
        this->EvalAsync(false);
        this->EvalFinished();

        const float *lower = this->lower.readOnlyHostPtr();
        const float *upper = this->upper.readOnlyHostPtr();
        const int *nbins = this->nbins.readOnlyHostPtr();
        const unsigned int *bins = this->bins->readOnlyHostPtr();
        const int *bin_stride = this->bin_stride.readOnlyHostPtr();
        const unsigned int norm = this->norm_buffer->readOnlyHostPtr()[this->norm_offset];

        TH1 *hist = 0;
        switch (this->nobservables) {
            case 1:
            hist = new TH1D("", "", nbins[0], lower[0], upper[0]);
            hist->Sumw2();
            for (int x=0; x < nbins[0]; x++) {
                int source_bin = x * bin_stride[0];
                int target_bin = hist->GetBin(x+1);
                hist->SetBinContent(target_bin, bins[source_bin] / this->bin_volume / norm);
                hist->SetBinError(target_bin, sqrt(bins[source_bin]) / this->bin_volume / norm);
            }
            break;

            case 2:
            hist = new TH2D("", "",
                            nbins[0], lower[0], upper[0],
                            nbins[1], lower[1], upper[1]);
            hist->Sumw2();
            for (int x=0; x < nbins[0]; x++) {
                for (int y=0; y < nbins[1]; y++) {
                    int source_bin = x * bin_stride[0] + y * bin_stride[1];
                    int target_bin = hist->GetBin(x+1, y+1);
                    hist->SetBinContent(target_bin, bins[source_bin] / this->bin_volume / norm);
                    hist->SetBinError(target_bin, sqrt(bins[source_bin]) / this->bin_volume / norm);
                }
            }
            break;

            case 3:
            hist = new TH3D("", "",
                            nbins[0], lower[0], upper[0],
                            nbins[1], lower[1], upper[1],
                            nbins[2], lower[2], upper[2]);
            hist->Sumw2();
            for (int x=0; x < nbins[0]; x++) {
                for (int y=0; y < nbins[1]; y++) {
                    for (int z=0; z < nbins[2]; z++) {
                        int source_bin = x * bin_stride[0] + y * bin_stride[1] + z * bin_stride[2];
                        int target_bin = hist->GetBin(x+1, y+1, z+1);
                        hist->SetBinContent(target_bin, bins[source_bin] / this->bin_volume / norm);
                        hist->SetBinError(target_bin, sqrt(bins[source_bin]) / this->bin_volume / norm);
                    }
                }
            }
            break;

            default:
                throw Error("Impossible EvalHist::CreateHistogram switch case!");
        }

        this->needs_optimization = orig_optimization_flag;

        return hist;
    }

    void EvalHist::Optimize()
    {
        OptimizeBin();
        OptimizeEval();
        needs_optimization = false;
    }

    void EvalHist::OptimizeBin()
    {
        // Only do this if running on GPU
        #ifdef __CUDACC__

        const int ngrid_sizes = 6;
        const int grid_sizes[ngrid_sizes] = { 2, 4, 8, 16, 32, 64 };
        const int nblock_sizes = 5;
        const int block_sizes[nblock_sizes] = { 32, 64, 128, 256, 512};
        const int nsamples = this->samples.size() / this->nfields;
        int nreps = 1;

        // Do more repetitions on small kernels to avoid being fooled by timing fluctuations
        if (nsamples < 1000)
            nreps = 1000;
        if (nsamples < 10000)
            nreps = 100;
        if (nsamples < 100000)
            nreps = 10;

        // Avoid picking large numbers of blocks due to timing fluctuations
        // by requiring the timing to be at least 10% better than the
        // current winner.  Since we try grid sizes in increasing order,
        // this favors smaller configurations, which will be better
        // for overlapping kernels.
        const double improvement_threshold = 0.9;

        // Handle missing systematics case
        int nsyst = 0;
        const SystematicDescriptor *syst_ptr = 0;
        if (this->syst) {
            nsyst = this->syst->size();
            syst_ptr = this->syst->readOnlyPtr();
        }

        // Force allocation of these buffers (since we are not calling the zero bin kernel first)
        this->bins->writeOnlyPtr(); 
        this->norm_buffer->writeOnlyPtr();

        // Benchmark all possible grid sizes
        float best_time = 1e9;
        TStopwatch timer;

        for (int igrid=0; igrid < ngrid_sizes; igrid++) {
            int grid_size = grid_sizes[igrid];

            for (int iblock=0; iblock < nblock_sizes; iblock++) {
                int block_size = block_sizes[iblock];

                // Skip obviously too-large launch configurations
                if (grid_size * block_size >= 2 * nsamples)
                    continue;

                timer.Start();
                for (int irep=0; irep < nreps; irep++) {
                    HEMI_KERNEL_LAUNCH(bin_samples, grid_size, block_size, 0, this->cuda_state->stream,
                      (int) this->samples.size(), this->samples.readOnlyPtr(), 
                      this->nobservables, this->nfields,
                      this->bin_stride.readOnlyPtr(), this->nbins.readOnlyPtr(),
                      this->lower.readOnlyPtr(), this->upper.readOnlyPtr(),
                      nsyst, syst_ptr,
                      this->param_buffer->readOnlyPtr() + this->param_offset, this->param_stride,
                      this->bins->ptr(), this->norm_buffer->writeOnlyPtr() + this->norm_offset);
                }
                checkCuda( cudaStreamSynchronize(this->cuda_state->stream) );
                timer.Stop();
                float this_time = timer.RealTime();

                //std::cerr << "Grid/Block tested:" << grid_size << "/" << block_size << "(" << this_time << " sec)\n";

                if (this_time < best_time * improvement_threshold) {
                    best_time = this_time;
                    this->bin_nblocks = grid_size;
                    this->bin_nthreads_per_block = block_size;
                }
            }
        }

        //#ifdef DEBUG
        std::cerr << "pdfz::EvalHist::OptimizeBin(): # of samples = " << nsamples
                  << " Grid/Block selected = " << this->bin_nblocks << "/" << this->bin_nthreads_per_block << "\n";
        //#endif
   
        #endif
    }


    void EvalHist::OptimizeEval()
    {
        // Only do this if running on GPU
        #ifdef __CUDACC__

        const int ngrid_sizes = 6;
        const int grid_sizes[ngrid_sizes] = { 2, 4, 8, 16, 32, 64 };
        const int nblock_sizes = 5;
        const int block_sizes[nblock_sizes] = { 32, 64, 128, 256, 512};
        const int npoints = this->read_bins->size();
        int nreps = 1;

        // Do more repetitions on small kernels to avoid being fooled by timing fluctuations
        if (npoints < 1000)
            nreps = 1000;
        if (npoints < 10000)
            nreps = 100;
        if (npoints < 100000)
            nreps = 10;

        // Avoid picking large numbers of blocks due to timing fluctuations
        // by requiring the timing to be at least 10% better than the
        // current winner.  Since we try grid sizes in increasing order,
        // this favors smaller configurations, which will be better
        // for overlapping kernels.
        const double improvement_threshold = 0.9;

        // Force allocation of these buffers (since we are not calling the zero bin kernel first)
        this->bins->writeOnlyPtr(); 
        this->norm_buffer->writeOnlyPtr();

        // Benchmark all possible grid sizes
        float best_time = 1e9;
        TStopwatch timer;

        for (int igrid=0; igrid < ngrid_sizes; igrid++) {
            int grid_size = grid_sizes[igrid];

            for (int iblock=0; iblock < nblock_sizes; iblock++) {
                int block_size = block_sizes[iblock];

                // Skip obviously too-large launch configurations
                if (grid_size * block_size >= 2 * npoints)
                    continue;

                timer.Start();
                for (int irep=0; irep < nreps; irep++) {
                    HEMI_KERNEL_LAUNCH(eval_pdf, grid_size, block_size, 0, this->cuda_state->stream,
                           (int) this->read_bins->size(), this->read_bins->readOnlyPtr(),
                           this->bins->readOnlyPtr(), this->norm_buffer->readOnlyPtr() + this->norm_offset,
                           this->bin_volume,
                           this->pdf_buffer->writeOnlyPtr() + this->pdf_offset, this->pdf_stride);
                }
                checkCuda( cudaStreamSynchronize(this->cuda_state->stream) );
                timer.Stop();
                float this_time = timer.RealTime();

                //std::cerr << "Grid/Block tested:" << grid_size << "/" << block_size << "(" << this_time << " sec)\n";

                if (this_time < best_time * improvement_threshold) {
                    best_time = this_time;
                    this->eval_nblocks = grid_size;
                    this->eval_nthreads_per_block = block_size;
                }
            }
        }

        //#ifdef DEBUG
        std::cerr << "pdfz::EvalHist::OptimizeEval(): # of points = " << npoints
                  << " Grid/Block selected = " << this->eval_nblocks << "/" << this->eval_nthreads_per_block << "\n";
        //#endif
   
        #endif
    }

    ///////////////////// EvalKernel ///////////////////////

} // namespace pdfz
