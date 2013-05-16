/**
 * \class NLL
 * \brief Evaluate the likelihood function
 *
 * This functor evaluates the negative log likelihood function. A lookup table
 * of the PDFs evaluated at each event (Pj(xi)) is constructed, and an
 * evaluation consists of looping over events and scaling the LUT entries by
 * the desired parameters.
 *
 * This time-consuming loop is run on a GPU using CUDA.
 */

#ifndef __NLL_H__
#define __NLL_H__

#include "signals.h"

class NLL {
  public:
    /**
     * Constructor
     *
     * \param signals set of Signals defining the PDFs
     * \param data TNtuple containing the experiment's dataset
     */
    NLL(const std::vector<Signal>& signals, TNtuple* data);

    virtual ~NLL();

    /**
     * Evaluate the NLL function at the given parameter values.
     *
     * \param _norms Normalizations at which to evaluate the function.
     *               These should be ordered the same as the signals.
     */
    double operator()(float* _norms);

  protected:
    /**
     * Build a lookup table of probabilities Pj(xi).
     *
     * \param signals set of Signals defining the PDFs
     * \param data TNtuple containing the experiment's dataset
     */
    static float* build_lut(const std::vector<Signal>& signals, TNtuple* data);

  private:
    size_t blocksize;  //!< number of threads per block
    size_t nblocks;  //!< number of blocks
    size_t nthreads;  //!< number of CUDA threads to run
    unsigned nsignals;  //!< number of signal dimensions
    unsigned nevents;  //!< number of events in dataset
    double* expectations;  //!< signals expectation values
    double* constraints;  //!< fractional constraints
    double* sums;  //!< output sums, host side
    double* sums_device;  //!< output sums, device side
    float* normalizations;  //!< normalizations, host side
    float* normalizations_device;  //!< normalizations, device side
    float* lut;  //!< lookup table, host side
    float* lut_device;  //!< lookup table, device side
};

#endif  // __NLL_H__

