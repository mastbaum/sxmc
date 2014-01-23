/**
 * \file mcmc.h
 *
 * Utilities for Markov Chain Monte Carlo distribution sampling.
 */

#ifndef __MCMC_H__
#define __MCMC_H__

#include <vector>
#include <cmath>
#include <string>
#include <cuda.h>
#include <hemi/hemi.h>

#include <sxmc/signals.h>
#include <sxmc/nll_kernels.h>
#include <sxmc/pdfz.h>

#ifdef __CUDACC__
#include <curand_kernel.h>
#endif

#ifndef __HEMI_ARRAY_H__
#define __HEMI_ARRAY_H__
#include <hemi/array.h>
#endif

class TNtuple;
class LikelihoodSpace;

/**
 * \class MCMC
 * \brief Markov Chain Monte Carlo simulator
 *
 * Given a set of signal PDFs and a dataset, random walk to map out the
 * likelihood space.
 */
class MCMC {
  public:
    /**
     * Constructor
     *
     * \param signals List of Signals defining the PDFs and expectations
     * \param systematics List of systematic parameter definitions
     * \param observables List of observables in the data
     */
    MCMC(const std::vector<Signal>& signals,
         const std::vector<Systematic>& systematics,
         const std::vector<Observable>& observables);

    /**
     * Destructor
     *
     * Free HEMI arrays
     */
    ~MCMC();

    /**
     * Perform walk.
     *
     * \param data TODO
     * \param weights TODO
     * \param nsteps Number of random-walk steps to take
     * \param burnin_fraction Fraction of initial steps to throw out
     * \param debug_mode If true, accept and save all steps
     * \param sync_interval How often to copy accepted from GPU to storage
     * \returns LikelihoodSpace built from samples
     */
    LikelihoodSpace* operator()(std::vector<float>& data,
                                std::vector<int>& weights,
                                unsigned nsteps,
                                float burnin_fraction,
                                const bool debug_mode=false,
                                unsigned sync_interval=10000);

  protected:
    /**
     * Evaluate the NLL function
     *
     * -logL = sum(Nj) + 1/2*sum((r-r')^2/s^2) - sum(log(sum(Nj*Pj(xi))))
     *
     * Nothing is returned -- the output stays in the nll array, so it can stay
     * on the device.
     *
     * This is done in three steps, to support running on the GPU:
     *
     *   1. Compute partial sums of chunks of events for the last term
     *   2. Total up partial sums from step 1
     *   3. Add normalization and other constraints with sum from step 2
     *
     * \param v Parameter vector at which to evaluate
     * \param nll Container for output NLL value
     * \param event_partial_sums Pre-allocated buffer for event term
     *                           calculation
     * \param event_total_sum Pre-allocated buffer for event term total
     */
    void nll(const float* lut, const int* dataweights, size_t nevents,
             const double* v, double* nll,
             double* event_partial_sums,
             double* event_total_sum);

  private:
    size_t nsignals;  //!< number of signal parameters
    size_t nsystematics;  //!< number of systematic parameters
    size_t nparameters;  //!< total number of parameters
    size_t nobservables;  //!< number of observables in data
    unsigned nnllblocks;  //!< number of cuda blocks for nll partial sums
    unsigned nllblocksize;  //!< size of cuda blocks for nll partial sums
    unsigned nnllthreads;  //!< number of threads for nll partial sums
    unsigned nreducethreads; //!< number of threads to use in partial sum
                             //!< reduction kernel
    unsigned blocksize;  //!< size of blocks for per-signal kernels
    unsigned nblocks;  //!< number of blocks for per-signal kernels
    std::string varlist;  //!< string identifier list for ntuple indexing
    hemi::Array<double>* parameter_means;  //!< parameter central values
    hemi::Array<double>* parameter_sigma;  //!< parameter Gaussian uncertainty
    hemi::Array<RNGState>* rngs;  //!< CURAND RNGs, ignored in CPU mode
    std::vector<std::string> parameter_names;  //!< string name of each param
    std::vector<pdfz::Eval*> pdfs;  //!< references to signal pdfs
};

#endif  // __MCMC_H__

