/** GPU/CPU-based PDF evaluation
 *
 * The subclasses of pdfz::Eval (EvalHist and EvalKernel) take an array of 
 * samples from an unknown multidimensional distribution and estimate the
 * probability density function over a finite volume.
 *
 * These classes are designed to be incorporated into a larger GPU-based
 * program using Hemi to manage GPU buffers.  As a result, the API is
 * designed to minimizes unnecessary copying of data between the CPU and GPU
 * by instead registering pointers to Hemi arrays before evaluation.
 * Additionally, to make better use of CUDA devices with compute 
 * capability 2.0 and greater, each evaluator object has its own CUDA
 * stream.  Kernels on different streams can run in parallel, so when
 * there are many PDFs to evaluate, the best approach is to call the 
 * EvalAsync() method on each evaluator object first, then call EvalFinished()
 * to block until the calculation is complete.
 *
 * The evaluator object computes the value of the PDF only at specific
 * "evaluation points".  As part of the evaluation step, the samples can be 
 * transformed on the fly to mimic various systematic uncertainties,
 * like scale, offset or resolution.  A parameter buffer is used to control
 * these systematic transformations.  In addition, a normalization buffer
 * is updated during evaluation to record the total number of samples
 * that were within the PDF domain after the systematic transformations 
 * were applied.  This can be used to calculate an efficiency correction
 * in certain kinds of applications.
 * 
 *
 * Example usage
 * -------------
 *
 *     const int nobs = 2;
 *     const int nfields = 3;
 *     const int nsamples = 100;
 *     const int neval_points = 20;
 *     const int nparams = 2; // two systematic parameters
 *
 *     std::vector<float> samples(nsamples * nfields);
 *     std::vector<float> lower(nobs);
 *     std::vector<float> upper(nobs);
 *     std::vector<float> eval_points(neval_points * nobs);
 *
 *     // Fill vectors with appropriate information here....
 *     // .........
 *   
 *     // CPU/GPU storage of evaluation inputs and outputs
 *     hemi::Array<float> parameters(nparams);
 *     hemi::Array<float> pdf_values(neval_points);
 *     hemi::Array<unsigned int> normalizations(1); // Only 1 PDF!
 *   
 *     // Settings for PDF evaluator
 *     std::vector<int> nbins(2, 10); // 10 bins in each dimension
 *   
 *     // Setup evaluator
 *     pdfz::EvalHist pdf(samples, nfields, nobs, lower, upper, nbins);
 *     pdf.SetEvalPoints(eval_points);
 *     pdf.SetPDFValueBuffer(&pdf_values, 0, 1);  // no offset, unit stride
 *     pdf.SetNormalizationBuffer(&normalizations, 0);
 *     pdf.SetParameterBuffer(parameters);
 *   
 *     // Add some systematics
 *     // Scaling systematic on observable 1 controlled by parameter 0
 *     pdf.AddSystematic(pdfz::ScaleSystematic(1, 0));
 *     // Resolution systematic on observable 0 (with true value in field 2)
 *     // controlled by parameter 1
 *     pdf.AddSystematic(pdfz::ResolutionSystmatic(0, 2, 1));
 *  
 *     // Do it!
 *     pdf.EvalAsync();
 *     pdf.EvalFinished();
 */
#ifndef __PDFZ_H__
#define __PDFZ_H__

#include <hemi/array.h>
#include <string>
#include <vector>

namespace pdfz {

    /**
    * \struct Error
    * \brief Generic PDF-related exception
    */
    struct Error
    {
        Error(const std::string &_msg) { msg = _msg; }

        /** Cause of exception */
        std::string msg;
    };


    /** 
    * \struct Systematic
    * \brief  Base struct used to describe a systematic uncertainty
    */
    struct Systematic
    {
        enum Type {
          SHIFT,
          SCALE,
          RESOLUTION_SCALE,
        };

        Type type;

        Systematic(Type _type) : type(_type) { }
        // Needed to make this a polymorphic type for dynamic_cast in Eval::AddSystmatic
        virtual ~Systematic() { }
    };


    /** 
     * \struct ShiftSystematic
     * \brief  Offset an observable.
     *
     * Transform x' = x + p
     */
    struct ShiftSystematic : public Systematic
    {
    // index of observable and systematic parameter
    ShiftSystematic(int _obs, int _par) :
        Systematic(SHIFT), obs(_obs), par(_par) { }

        int obs, par;
    };


    /**
     * \struct ScaleSystematic
     * \brief Rescale an observable.
     *
     * Transform x' = x * (1 + p)
     */
    struct ScaleSystematic : public Systematic
    {
    // index of observable and systematic parameter
    ScaleSystematic(int _obs, int _par) :
        Systematic(SCALE), obs(_obs), par(_par) { }

        int obs, par;
    };


    /**
     * \struct ResolutionScaleSystematic
     * \brief  Fractionally alter the resolution of an observable by rescaling its distance
     *         from a "true" value.
     *
     * Transform x' = (1 + p) * (x - x_true) + x_true
     */
    struct ResolutionScaleSystematic : public Systematic
    {
    // index of observable, "true" observable, and systematic parameter
    ResolutionScaleSystematic(int _obs, int _true_obs, int _par) :
        Systematic(RESOLUTION_SCALE), obs(_obs), true_obs(_true_obs), par(_par) { }

        int obs, true_obs, par;
    };


    struct CudaState; // Hide CUDA-related state from callers
    struct SystematicDescriptor; // Also hide representation of systematics

    class Eval
    {
    public:
        /** Create a PDF evaluator for the PDF samples stored in ``samples``.
            
            ``samples`` is assumed to be a flattened 2D array in row-major order 
            (x1,y1,x2,y2,etc), where each row is ``nfields`` in length and the
            first ``nobservables`` elements of each row represent the dimensions
            of the PDF domain.  Extra fields, if any, may be used to control
            systematics.  The ``lower`` and ``upper`` arrays (size()==nobs) set 
            the lower and upper bounds for each observable.

            Raises pdfz::Error if samples.size() is not divisible by nfields, if
            nobservables > nfields, or if nobservables == 0. 
        */
        Eval(const std::vector<float> &samples, int nfields, int nobservables,
             const std::vector<float> &lower, const std::vector<float> &upper);

        virtual ~Eval();


        /** Set the points where the PDF will be evaluated.

            ``points`` is a flattened 2D array in row-major order
            (x1,y1,x2,y2,etc), where each row is this->GetNobs() in length.
            This method does non-trivial calculation, so it should not be
            called in performance-critical loops.
        */
        virtual void SetEvalPoints(const std::vector<float> &points)=0;


        /** Set the output array where the PDF values will be written for each point.

            To allow PDF values to be written into a larger array, this method
            saves a pointer to a hemi:Array to be used during evaluation.  The
            PDF values will be normalized such that the integral of the PDF
            over the range given by lower and upper in the constructor is 1.

            At evaulation time, the PDF evaluated at t_i will be written to:
                output[offset + i * stride]
        */
        virtual void SetPDFValueBuffer(hemi::Array<float> *output, int offset=0, int stride=1);


        /** Set the output array where the PDF normalization will be written for each point.

            To allow normalization values to be written into a larger array, this method
            saves a pointer to a hemi:Array to be used during evaluation.  The 
            units of the normalization is "number of events", so if all samples are
            fully contained within the upper and lower bounds of the PDF after systematic
            transformations, the normalization factor written during evaluation will
            be equal to samples.size().

            At evaulation time, the normalization will be written to:
                norm[offset]
        */
        virtual void SetNormalizationBuffer(hemi::Array<unsigned int> *norm, int offset=0);


        /** Set the storage buffer that will be read for systematic parameters.

            At evaluation time, the systematic parameter j will be read from:
                params[offset + j * stride]
        */
        virtual void SetParameterBuffer(hemi::Array<float> *params, int offset=0, int stride=1);


        /** Add a systematic transformation to this PDF */
        virtual void AddSystematic(const Systematic &syst);


        /** Launch evaluation of the PDF at all the points given in the last call to
            SetEvalPoints() using the systematic parameters read from the
            parameter buffer specified in SetParameterBuffer().

            PDF values are written to the output array given in SetOutputBuffer().

            Note that on systems with GPUs, this function returns before the evaluation
            completes.  The output and normalization arrays should not be read, and the 
            parameter buffer should not be written to until EvalFinished() is called.
        */
        virtual void EvalAsync()=0;

        /** Wait until the evaluation launched in EvalAsync() has completed.

        Once this function returns, the output and normalization arrays will be set 
        and can be read, and the parameter buffer can be written to without creating
        a race condition.
        */
        virtual void EvalFinished()=0;

    protected:
        int nfields;
        int nobservables;

        hemi::Array<float> lower;
        hemi::Array<float> upper;

        hemi::Array<float> *pdf_buffer;
        int pdf_offset;
        int pdf_stride;

        hemi::Array<unsigned int> *norm_buffer;
        int norm_offset;

        hemi::Array<float> *param_buffer;
        int param_offset;
        int param_stride;

        hemi::Array<SystematicDescriptor> *syst;

        CudaState *cuda_state;
    };


    /** Evaluate a PDF using an N-dimensional histogram */
    class EvalHist : public Eval
    {
    public:

        /** See Eval::Eval() for description of parameters.  
            ``nbins`` gives the number of histogram bins for each
            observable dimension.

            Raises all exceptions of Eval(), as well as pdfz::Error if 
            nbins.size() != nobs.

            If optimize is set to true (default), then the CUDA block
            configuration will be optimized the first time EvalAsync() is called.
        */
        EvalHist(const std::vector<float> &samples, int nfields, int nobservables,
                 const std::vector<float> &lower, const std::vector<float> &upper,
                 const std::vector<int> &nbins, bool optimize=true);

        virtual ~EvalHist();
        virtual void SetEvalPoints(const std::vector<float> &points);

        /** Brute force tests a bunch of CUDA configurations to find the best one */
        virtual void Optimize();

        virtual void EvalAsync();
        virtual void EvalFinished();

    protected:
        hemi::Array<float> samples;
        hemi::Array<float> *eval_points;
        hemi::Array<int> nbins;
        hemi::Array<int> bin_stride;
        int total_nbins;
        hemi::Array<unsigned int> *bins;

        int nthreads_per_block;
        int nblocks;

        bool needs_optimization;
    };


    /** Evaluate a PDF using a kernel density estimator */
    class EvalKernel : public Eval
    {
        /** See Eval::Eval() for description of parameters.  
            ``bandwidth scale`` gives rescales the bandwidths in each dimension
            by the given value.  Pass an array of 1.0f to use the default
            bandwidth calculation.

            Raises all exceptions of PDFEval(), as well as pdfz::Error if 
            bandwidth_scale.size() != nobs.
        */
        EvalKernel(const std::vector<float> &samples, int nfields, int nobservables,
                   const std::vector<float> &lower, const std::vector<float> &upper,
                   const std::vector<float> &bandwidth_scale);

        virtual ~EvalKernel();
        virtual void SetEvalPoints(const std::vector<float> &points);
        virtual void EvalAsync();
        virtual void EvalFinished();
    };

} // end namespace pdfz

#endif // __PDFZ_H__
