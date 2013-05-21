#ifndef __MCMC_H__
#define __MCMC_H__

#include <vector>
#include <cmath>
#include <cuda.h>
#include <hemi/hemi.h>
#include "signals.h"

#ifdef __CUDACC__
#include <curand_kernel.h>
#endif

#ifndef __HEMI_ARRAY_H__
#define __HEMI_ARRAY_H__
#include <hemi/array.h>
#endif

#ifdef __CUDACC__
typedef curandState RNGState;
#else
typedef int RNGState;  // ignored by CPU code
#endif

class TNtuple;

#ifdef __CUDACC__
/**
 * Initialize device-side RNGs.
 *
 * Generators all have the same seed but a different offset in the sequence.
 *
 * \param nthreads Number of threads (same as the number of states)
 * \param seed Random seed shared by all generators
 * \param state Array of CUDA RNG states
 */
__global__ void init_device_rngs(int nthreads, unsigned long long seed,
                                 curandState* state);
#endif


/**
 * Pick a new position distributed around the given one.
 *
 * Uses CURAND XORWOW generator on GPU, or ROOT's gRandom on the CPU.
 *
 * \param nthreads Number of threads == length of vectors
 * \param rng CUDA RNG states, ignored on CPU
 * \param sigma Standard deviation to sample
 * \param current_vector Vector of current parameters
 * \param proposed_vector Output vector of proposed parameters
 */
HEMI_KERNEL(pick_new_vector)(int nthreads, RNGState* rng, float sigma,
                             float* current_vector, float* proposed_vector);


/**
 * NLL Part 1
 *
 * Calculate -sum(log(sum(Nj * Pj(xi)))) contribution to NLL.
 *
 * \param lut Pj(xi) lookup table
 * \param pars Event rates (normalizations) for each signal
 * \param ne Number of events in the data
 * \param ns Number of signals
 * \param sums Output sums for subsets of events
 */
HEMI_KERNEL(nll_event_chunks)(const float* lut, const float* pars,
                              const size_t ne, const size_t ns,
                              double* sums);


/**
 * NLL Part 2
 *
 * Total up the partial sums from Part 1
 *
 * \param
 */
HEMI_KERNEL(nll_event_reduce)(const size_t nthreads, const double* sums,
                              double* total_sum);


/**
 * NLL Part 3
 *
 * Calculate overall normalization and constraints contributions to NLL, add
 * in the event term to get the total.
 *
 * \param ns Number of signals
 * \param pars Event rates (normalizations) for each signal
 * \param expectations Expected rates for each signal
 * \param constraints Fractional constraints for each signal
 * \param events_total Sum of event term contribution
 * \param nll The total NLL
 */
HEMI_KERNEL(nll_total)(const size_t ns, const float* pars,
                       const float* expectations,
                       const float* constraints,
                       const double* events_total,
                       float* nll);


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
    MCMC(std::vector<Signal> signals, TNtuple* data);

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
    void nll(hemi::Array<float>& v, hemi::Array<float>& nll,
             hemi::Array<double>& event_partial_sums,
             hemi::Array<double>& event_total_sum);

  private:
    unsigned nsignals;
    unsigned nevents;
    unsigned nnllblocks;
    unsigned nllblocksize;
    unsigned nnllthreads;
    unsigned blocksize;
    unsigned nblocks;
    hemi::Array<float>* expectations;
    hemi::Array<float>* constraints;
    hemi::Array<float>* lut;  // Event/PDF probability lookup table
    hemi::Array<RNGState>* rngs;  // CURAND RNGs, ignored in CPU mode
};

#endif  // __MCMC_H__

