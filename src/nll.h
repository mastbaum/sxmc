#ifndef __NLL_H__
#define __NLL_H__

#include "signals.h"

// functor to evaluate the negative log likelihood on the gpu
class NLL {
  public:
    NLL(const std::vector<Signal>& signals, TNtuple* data);
    virtual ~NLL();

    // evaluate the NLL function at the given parameter values
    double operator()(float* _norms);

  protected:
    // build a lookup table of probabilities Pj(xi)
    static float* build_lut(const std::vector<Signal>& signals, TNtuple* data);

  private:
    size_t blocksize;  // number of threads per block
    size_t nblocks;  // number of blocks
    size_t nthreads;  // number of CUDA threads to run
    unsigned nsignals;  // number of signal dimensions
    unsigned nevents;  // number of events in dataset
    double* expectations;  // signals expectation values
    double* constraints;  // fractional constraints
    double* sums;  // output sums, host side
    double* sums_device;  // output sums, device side
    float* normalizations;  // normalizations, host side
    float* normalizations_device;  // normalizations, device side
    float* lut;  // lookup table, host side
    float* lut_device;  // lookup table, device side
};

#endif  // __NLL_H__

