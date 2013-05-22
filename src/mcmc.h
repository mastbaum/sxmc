#ifndef __MCMC_H__
#define __MCMC_H__

#include <vector>
#include <cmath>
#include <string>
#include <cuda.h>
#include <hemi/hemi.h>
#include "signals.h"
#include "nll_kernels.h"

#ifdef __CUDACC__
#include <curand_kernel.h>
#endif

#ifndef __HEMI_ARRAY_H__
#define __HEMI_ARRAY_H__
#include <hemi/array.h>
#endif

class TNtuple;

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
     * \param data The dataset as a TNtuple with fields "e:r"
     */
    MCMC(const std::vector<Signal>& signals, TNtuple* data);

    /**
     * Destructor
     *
     * Free HEMI arrays
     */
     ~MCMC();

    /**
     * Perform walk.
     *
     * \param nsteps Number of random-walk steps to take
     * \param burnin_fraction Fraction of initial steps to throw out
     * \param sync_interval How often to copy accepted from GPU to storage
     * \returns TNtuple containing accepted points and their likelihoods
     */
    TNtuple* operator()(unsigned nsteps, float burnin_fraction,
                        unsigned sync_interval=1000);

    /**
     * Build Pj(xi) lookup table
     *
     * Build a table with nevents rows and nsignals columns, containing the
     * value of normalized PDF j evaluted at event i.
     *
     * \param signals Signals defining the PDFs
     * \param data TNtuple defining the dataset
     * \returns HEMI array of probabilities, indexed as noted above
     */
    hemi::Array<float>* build_lut(const std::vector<Signal>& signals,
                                        TNtuple* data);

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
     * \param event_partial_sums Pre-allocated buffer for event term calculation
     * \param event_total_sum Pre-allocated buffer for event term total
     */
    void nll(const float* v, double* nll,
             double* event_partial_sums,
             double* event_total_sum);

  private:
    unsigned nsignals;  //!< number of signals
    unsigned nevents;  //!< number of events
    unsigned nnllblocks;  //!< number of cuda blocks for nll partial sums
    unsigned nllblocksize;  //!< size of cuda blocks for nll partial sums
    unsigned nnllthreads;  //!< number of threads for nll partial sums
    unsigned nreducethreads; //!< number of threads to use in partial sum reduction kernel
    unsigned blocksize;  //!< size of blocks for per-signal kernels
    unsigned nblocks;  //!< number of blocks for per-signal kernels
    std::string varlist;  //!< string identifier list for ntuple indexing
    hemi::Array<float>* expectations;  //!< signal rate expectation values
    hemi::Array<float>* constraints;  //!< signal rate gaussian constraints
    hemi::Array<float>* lut;  //!< Event/PDF probability lookup table
    hemi::Array<RNGState>* rngs;  //!< CURAND RNGs, ignored in CPU mode
    std::vector<std::string> signal_names;  //!< string name of each signal
};

#endif  // __MCMC_H__

