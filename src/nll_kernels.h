#ifndef __NLL_KERNELS_H__
#define __NLL_KERNELS_H__

/**
 * \file nll_kernels.h
 * \brief CUDA/HEMI kernels supporting NLL calculation
*/

#include <cuda.h>
#include <hemi/hemi.h>

#ifdef __CUDACC__
#include <curand_kernel.h>
#endif

#ifndef __HEMI_ARRAY_H__
#define __HEMI_ARRAY_H__
#include <hemi/array.h>
#endif

/**
 * \typedef RNGState
 * \brief Defines RNG for CURAND, ignored in CPU mode
*/
#ifdef __CUDACC__
typedef curandStateXORWOW RNGState;
#else
typedef int RNGState;
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
 * \param jump_width Standard deviations to sample for each dimension
 * \param current_vector Vector of current parameters
 * \param proposed_vector Output vector of proposed parameters
*/
HEMI_KERNEL(pick_new_vector)(int nthreads, RNGState* rng,
                             const float* jump_width,
                             const double* current_vector,
                             double* proposed_vector);


/**
 * Decide whether to accept a random MCMC step.
 *
 * Compare likelihoods of current and proposed parameter vectors. Store each
 * step in a buffer which can be flushed periodically, minimizing transfer
 * overhead.
 *
 * The step buffer is an (Nsignals + 1 x Nsteps) matrix, where the last column
 * contains the likelihood value.
 *
 * \param rng Random-number generator states, used in GPU mode only
 * \param nll_current The NLL of the current parameters
 * \param nll_proposed the NLL of the proposed parameters
 * \param v_current The current parameters
 * \param v_proposed The proposed parameters
 * \param nparameters The number of parameters
 * \param accepted Number of accepted steps
 * \param counter The number of steps in the buffer
 * \param jump_buffer The step buffer
*/
HEMI_KERNEL(jump_decider)(RNGState* rng, double* nll_current,
                          const double* nll_proposed, double* v_current,
                          const double* v_proposed, unsigned nparameters,
                          int* accepted, int* counter, float* jump_buffer);


/**
 * NLL Part 1
 *
 * Calculate -sum(log(sum(Nj * Pj(xi)))) contribution to NLL.
 *
 * \param lut Pj(xi) lookup table
 * \param dataweights TODO
 * \param pars Event rates (normalizations) for each signal
 * \param ne Number of events in the data
 * \param ns Number of signals
 * \param norms The number of events in the PDF range (with systematics)
 * \param norms_nominal The number of events nominally in the PDF
 * \param sums Output sums for subsets of events
*/
HEMI_KERNEL(nll_event_chunks)(const float* lut, const int* dataweights,
                              const double* pars, const size_t ne,
                              const size_t ns,
                              const unsigned* norms,
                              const unsigned* norms_nominal,
                              double* sums);


/**
 * NLL Part 2
 *
 * Total up the partial sums from Part 1
 *
 * \param nthreads Number of threads == number of sums to total
 * \param sums The partial sums
 * \param total_sum Output: the total sum
*/
HEMI_KERNEL(nll_event_reduce)(const size_t nthreads, const double* sums,
                              double* total_sum);


/**
 * NLL Part 3
 *
 * Calculate overall normalization and constraints contributions to NLL, add
 * in the event term to get the total.
 *
 * \param nparameters The number of parameters
 * \param pars Parameters, normalizations then systematics
 * \param nsignals Number of signal parameters
 * \param means Expected rates and means of systematics
 * \param sigmas Gaussian constraint sigma, same units as means
 * \param events_total Sum of event term contribution
 * \param norms The number of events in the PDF range (with systematics)
 * \param norms_nominal The number of events nominally in the PDF
 * \param nll The total NLL
*/
HEMI_KERNEL(nll_total)(const size_t nparameters, const double* pars,
                       const size_t nsignals,
                       const double* means,
                       const double* sigmas,
                       const double* events_total,
                       const unsigned* norms,
                       const unsigned* norms_nominal,
                       double* nll);


/**
 * All-in-one MCMC step kernel.
 *
 * Combines the operations of computing the final NLL, picking a new parameter
 * vector, and possibly jumping, all in one kernel to reduce launch overhead.
 *
 * \param npartial_sums The number of partial sums of event terms to add up
 * \param sums Partial sums from event terms
 * \param ns The number of signals
 * \param means Expected rates and means of systematics
 * \param sigmas Gaussian constraint sigma, same units as means
 * \param rng Random-number generators
 * \param nll_current The NLL at the current step
 * \param nll_proposed The NLL at the proposed step
 * \param v_current The current parameter vector
 * \param v_proposed The proposed parameter vector
 * \param accepted The number of accepted steps
 * \param counter The number of steps in the jump buffer
 * \param jump_buffer The buffer of steps (vectors and likelihoods)
 * \param nparameters The number of parameters (dimensions in the L space)
 * \param jump_width The jump distribution widths in each dimension
 * \param norms The number of events in the PDF range (with systematics)
 * \param norms_nominal The number of events nominally in the PDF
 * \param debug_mode Enable debugging mode, where every step is accepted
*/
HEMI_KERNEL(finish_nll_jump_pick_combo)(const size_t npartial_sums,
                                        const double* sums, const size_t ns,
                                        const double* means,
                                        const double* sigmas,
                                        RNGState* rng,
                                        double *nll_current,
                                        double *nll_proposed,
                                        double *v_current, double *v_proposed,
                                        int* accepted, int* counter,
                                        float* jump_buffer, int nparameters,
                                        const float* jump_width,
                                        const unsigned* norms,
                                        const unsigned* norms_nominal,
                                        const bool debug_mode=false);

#endif  // __NLL_KERNELS_H__

