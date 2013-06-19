#include "pdfz.h"
#include <math.h>
#include <cuda.h>
#include <math_constants.h> // CUDA header

namespace pdfz {

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

    // If not compiling with CUDA, alias the CUDA stream type to something harmless.
    #ifndef __CUDACC__
    typedef int cudaStream_t;
    #endif


    // Hidden CUDA state for evaluator.  Not visible to Eval class users.
    struct CudaState
    {
        cudaStream_t stream;
    };


    //HEMI_ALIGN(8)
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
                       const std::vector<int> &_nbins) :
        Eval(_samples, nfields, nobservables, lower, upper),
        samples(_samples.size(), false), eval_points(0), 
        nbins(_nbins.size(), true), bin_stride(_nbins.size(), true), bins(0)
    {
        if ( (int) _nbins.size() != nobservables)
            throw Error("Size of nbins array must be same as number of observables.");

        this->samples.copyFromHost(&_samples.front(), _samples.size());
        this->nbins.copyFromHost(&_nbins.front(), _nbins.size());

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


        this->nthreads_per_block = 256;
        if (_samples.size() > 100000)
            this->nblocks = 64;
        else if (_samples.size() > 10000)
            this->nblocks = 16;
        else
            this->nblocks = 4;
    }

    EvalHist::~EvalHist()
    {
        delete this->eval_points;
        delete this->bins;
    }

    void EvalHist::SetEvalPoints(const std::vector<float> &points)
    {
        if (points.size() % this->nobservables != 0)
            throw Error("Number of entries in evaluation points array not divisible by number of observables.");

        delete this->eval_points;
        this->eval_points = new hemi::Array<float>(points.size(), false);
        this->eval_points->copyFromHost(&points.front(), points.size());
    }

    ///// EvalHist kernels
    HEMI_DEV_CALLABLE_INLINE
    void apply_systematic(const SystematicDescriptor *syst, float *fields, const float *parameters)
    {
        switch (syst->type) {
            case Systematic::SHIFT:
            fields[syst->obs] += parameters[syst->par];
            break;
            case Systematic::SCALE:
            fields[syst->obs] *= (1 + parameters[syst->par]);
            break;
            case Systematic::RESOLUTION_SCALE:
            fields[syst->obs] *= (1 + parameters[syst->par]);            
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

    HEMI_KERNEL(bin_samples)(int nsamples, const float *samples,
                             const int nobs, const int nfields, 
                             const int *bin_stride, const int *nbins,
                             const float *lower, const float *upper,
                             const int nsyst, const SystematicDescriptor *syst,
                             const float *parameters,
                             unsigned int *bins, unsigned int *norm)
    {
        int offset = hemiGetElementOffset();
        int stride = hemiGetElementStride();
        float *field_buffer = new float[nfields];

        for (int isample=offset; isample < nsamples; isample += stride) {
            bool in_pdf_domain = true;
            int bin_id = 0;

            // Copy fields
            for (int ifield=0; ifield < nfields; ifield++)
                field_buffer[ifield] = samples[isample * nfields + ifield];

            // Apply systematics
            for (int isyst=0; isyst < nsyst; isyst++)
                apply_systematic(syst + isyst, field_buffer, parameters);

            // Compute histogram bin
            for (int iobs=0; iobs < nobs; iobs++) {
                float element = field_buffer[iobs];
                // Throw out this event if outside of PDF domain
                if (element < lower[iobs] || element >= upper[iobs]) {
                    in_pdf_domain = false;
                    break;
                }

                bin_id += (int)( ((element - lower[iobs]) / (upper[iobs] - lower[iobs])) * nbins[iobs] ) * bin_stride[iobs];
            }

            // Add to histogram if sample in PDF domain
            if (in_pdf_domain) {
                #ifdef HEMI_DEV_CODE
                atomicAdd(bins + bin_id, 1);
                atomicAdd(norm, 1);
                #else
                bins[bin_id] += 1;
                *norm += 1;
                #endif
            }
        }

        delete[] field_buffer;
    }

    HEMI_KERNEL(eval_pdf)(int npoints, const float *points, int point_stride,
                          int ndims, const int *bin_stride, const int *nbins,
                          const float *lower, const float *upper,
                          const unsigned int *bins, const unsigned int *norm,
                          float *output, int output_stride)
    {
        int offset = hemiGetElementOffset();
        int stride = hemiGetElementStride();

        for (int ipoint=offset; ipoint < npoints; ipoint += stride) {
            bool in_pdf_domain = true;
            int bin_id = 0;
            float bin_volume = 1.0f;

            for (int idim=0; idim < ndims; idim++) {
                int ielement = point_stride * ipoint + idim;
                float element = points[ielement];

                // Throw out this event if outside of PDF domain
                if (element < lower[idim] || element >= upper[idim]) {
                    in_pdf_domain = false;
                    break;
                }

                float span = upper[idim] - lower[idim];
                bin_volume *= span / nbins[idim];
                bin_id += (int)( ((element - lower[idim]) / span) * nbins[idim] ) * bin_stride[idim];
            }

            float pdf_value = 0.0f;
            if (in_pdf_domain) {
                pdf_value = bins[bin_id] / bin_volume / *norm;
            } else {
                pdf_value = nanf(0);
            }
            output[output_stride * ipoint] = pdf_value;
        }
    }
    ///// End EvalHist kernels

    void EvalHist::EvalAsync()
    {
        // Handle no systematics case
        int nsyst = 0;
        const SystematicDescriptor *syst_ptr = 0;
        if (this->syst) {
            nsyst = this->syst->size();
            syst_ptr = this->syst->readOnlyPtr();
        }

        HEMI_KERNEL_LAUNCH(zero_hist, this->nblocks, this->nthreads_per_block, 0, this->cuda_state->stream,
                           this->total_nbins, this->bins->writeOnlyPtr(), this->norm_buffer->writeOnlyPtr() + this->norm_offset);
        HEMI_KERNEL_LAUNCH(bin_samples, this->nblocks, this->nthreads_per_block, 0, this->cuda_state->stream,
                           (int) this->samples.size(), this->samples.readOnlyPtr(), 
                           this->nobservables, this->nfields,
                           this->bin_stride.readOnlyPtr(), this->nbins.readOnlyPtr(),
                           this->lower.readOnlyPtr(), this->upper.readOnlyPtr(),
                           nsyst, syst_ptr,
                           this->param_buffer->readOnlyPtr(),
                           this->bins->ptr(), this->norm_buffer->writeOnlyPtr() + this->norm_offset);
        HEMI_KERNEL_LAUNCH(eval_pdf, this->nblocks, this->nthreads_per_block, 0, this->cuda_state->stream,
                           (int) this->eval_points->size(), this->eval_points->readOnlyPtr(), this->nobservables,
                           this->nobservables, this->bin_stride.readOnlyPtr(), this->nbins.readOnlyPtr(),
                           this->lower.readOnlyPtr(), this->upper.readOnlyPtr(), 
                           this->bins->readOnlyPtr(), this->norm_buffer->readOnlyPtr() + this->norm_offset,
                           this->pdf_buffer->writeOnlyPtr() + this->pdf_offset, this->pdf_stride);
    }

    void EvalHist::EvalFinished()
    {
        #ifdef __CUDACC__
        checkCuda( cudaStreamSynchronize(this->cuda_state->stream) );
        #endif
    }


    ///////////////////// EvalKernel ///////////////////////

} // namespace pdfs