#include "pdfz.h"
#include <math.h>
#include <cuda.h>
#include <math_constants.h> // CUDA header

namespace pdfz {

    // If not compiling with CUDA, alias the CUDA stream type to something harmless.
    #ifndef __CUDACC__
    typedef int cudaStream_t;
    #endif


    // Hidden CUDA state for evaluator.  Not visible to Eval class users.
    struct CudaState
    {
        cudaStream_t stream;
    };

    Eval::Eval(const std::vector<float> &_samples, int _nfields, int _nobservables,
               const std::vector<float> &_lower, const std::vector<float> &_upper) :
        nfields(_nfields), nobservables(_nobservables), 
        lower(_lower.size(), true), upper(_upper.size(), true)
    {
        if (_samples.size() % _nfields != 0)
            throw Error("Length of samples array is not divisible by number of fields.");

        if (_nobservables == 0 )
            throw Error("Number of observables in PDF is zero.");

        if (_nobservables > nfields )
            throw Error("Number of observables cannot be greater than number of fields.");

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
    }

    void Eval::SetPDFBuffer(hemi::Array<float> *output, int offset, int stride)
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

    }


    ///////////////////// EvalHist ///////////////////////

    EvalHist::EvalHist(const std::vector<float> &_samples, int nfields, int nobservables,
                       const std::vector<float> &lower, const std::vector<float> &upper,
                       const std::vector<int> &_nbins) :
        Eval(_samples, nfields, nobservables, lower, upper),
        samples(_samples.size(), false), eval_points(0), 
        nbins(_nbins.size(), true), bin_stride(_nbins.size(), true), bins(0)
    {
        if (_nbins.size() == 0)
            throw Error("Cannot create EvalHist object with 0 dimensions");

        this->samples.copyFromHost(&_samples.front(), _samples.size());
        this->nbins.copyFromHost(&_nbins.front(), _nbins.size());

        // Compute stride for row-major order storage
        const int ndims = (int) this->bin_stride.size();
        int *bin_stride_host = this->bin_stride.hostPtr();

        bin_stride_host[ndims - 1] = 1;
        for (int i=ndims-2; i >= 0; i--)
            bin_stride_host[i] = _nbins[i+1] * bin_stride_host[i+1];

        this->total_nbins = bin_stride_host[0] * _nbins[0];
        this->bins = new hemi::Array<unsigned int>(this->total_nbins, true);


        this->nthreads_per_block = 256;
        this->nblocks = 64;
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
    HEMI_KERNEL(zero_hist)(int total_nbins, unsigned int *bins, unsigned int *norm)
    {
        int offset = hemiGetElementOffset();
        int stride = hemiGetElementStride();

        if (offset == 0)
            *norm = 0;

        for (int i=offset; i < total_nbins; i += stride)
            bins[i] = 0;
    }

    HEMI_KERNEL(bin_samples)(int nsamples, const float *samples, int sample_stride,
                             int ndims, const int *bin_stride, const int *nbins,
                             const float *lower, const float *upper,
                             unsigned int *bins, unsigned int *norm)
    {
        int offset = hemiGetElementOffset();
        int stride = hemiGetElementStride();

        for (int isample=offset; isample < nsamples; isample += stride) {
            bool in_pdf_domain = true;
            int bin_id = 0;

            for (int idim=0; idim < ndims; idim++) {
                int ielement = sample_stride * isample + idim;
                float element = samples[ielement];

                // Throw out this event if outside of PDF domain
                if (element < lower[idim] || element >= upper[idim]) {
                    in_pdf_domain = false;
                    break;
                }

                bin_id += (int)( ((element - lower[idim]) / (upper[idim] - lower[idim])) * nbins[idim] ) * bin_stride[idim];
            }

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
                bin_volume *= span;
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
        HEMI_KERNEL_LAUNCH(zero_hist, this->nblocks, this->nthreads_per_block, 0, this->cuda_state->stream,
                           this->total_nbins, this->bins->writeOnlyPtr(), this->norm_buffer->writeOnlyPtr() + this->norm_offset);
        HEMI_KERNEL_LAUNCH(bin_samples, this->nblocks, this->nthreads_per_block, 0, this->cuda_state->stream,
                           (int) this->samples.size(), this->samples.readOnlyPtr(), this->nfields,
                           this->nobservables, this->bin_stride.readOnlyPtr(), this->nbins.readOnlyPtr(),
                           this->lower.readOnlyPtr(), this->upper.readOnlyPtr(), 
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