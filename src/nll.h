#ifndef __NLL_H__
#define __NLL_H__

#include "signals.h"

// functor to evaluate the negative log likelihood on the gpu
class NLL {
  public:
    NLL(const std::vector<Signal>& signals, TNtuple* data);
    virtual ~NLL();

    // evaluate the NLL function at the given parameter values
    double operator()(double* _norms);

  protected:
    // build a lookup table of probabilities Pj(xi)
    static float* build_lut(const std::vector<Signal>& signals, TNtuple* data);

  private:
    unsigned nsignals;  // number of signal dimensions
    unsigned nevents;  // number of events in dataset
    double* expectations;  // signals expectation values
    double* constraints;  // fractional constraints
    float* lut;  // lookup table, host side
    float* lut_device;  // lookup table, device side
};

#endif  // __NLL_H__

