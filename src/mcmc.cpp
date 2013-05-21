#include <iostream>
#include <vector>
#include <cmath>
#include <hemi/hemi.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TRandom.h>
#include "mcmc.h"
#include "signals.h"

#ifdef __CUDACC__
#include <curand_kernel.h>
#endif

#ifndef __HEMI_ARRAY_H__
#define __HEMI_ARRAY_H__
#include <hemi/array.h>
#endif

#ifdef __CUDACC__
__global__ void init_device_rngs(int nthreads, unsigned long long seed,
                                 curandState* state) {
  int idx = hemiGetElementOffset();
  if (idx > nthreads) {
    return;
  }
  curand_init(seed, idx, 0, &state[idx]);
}
#endif

HEMI_KERNEL(pick_new_vector)(int nthreads, curandState* rng, float sigma,
                             const float* current_vector,
                             float* proposed_vector) {
  int idx = hemiGetElementOffset();
  if (idx > nthreads) {
    return;
  }
#ifdef HEMI_DEV_CODE
  float u = curand_normal(&rng[idx]);
  proposed_vector[idx] = current_vector[idx] + sigma * u;
#else
  proposed_vector[idx] = gRandom->Gaus(current_vector[idx], sigma);
#endif
}


HEMI_KERNEL(nll_event_chunks)(const float* lut, const float* pars,
                              const size_t ne, const size_t ns,
                              double* sums) {
  int offset = hemiGetElementOffset();
  int stride = hemiGetElementStride();

  double sum = 0;
  for (int i=offset; i<(int)ne; i+=stride) {
    double s = 0;
    for (size_t j=0; j<ns; j++) {
      s += pars[j] * lut[i * ns + j];
    }
    sum += log(s);
  }

  sums[offset] = sum;
}


HEMI_KERNEL(nll_event_reduce)(const size_t nthreads, const double* sums,
                             double* total_sum) {
  float sum = 0;
  for (size_t i=0; i<nthreads; i++) {
    sum += sums[i];
  }
  *total_sum = sum;
}


HEMI_KERNEL(nll_total)(const size_t ns, const float* pars,
                       const float* expectations,
                       const float* constraints,
                       const double* events_total,
                       float* nll) {
  int idx = hemiGetElementOffset();

  // total from sum over events, once
  if (idx == 0) {
#ifdef HEMI_DEV_CODE
    atomicAdd(nll, *events_total);
#else
    *nll = *nll + *events_total;
#endif
  }

  // normalization constraints
#ifdef HEMI_DEV_CODE
  atomicAdd(nll, pars[idx]);
#else
  *nll = *nll + pars[idx];
#endif

  // gaussian constraints
  if (constraints[idx] != 0) {
    float x = (pars[idx] / expectations[idx] - 1.0) / constraints[idx];
#ifdef HEMI_DEV_CODE
    atomicAdd(nll, x * x);
#else
    *nll = *nll + x * x;
#endif
  }
}


MCMC::MCMC(std::vector<Signal> signals, TNtuple* data) {
  this->nsignals = signals.size();
  this->nevents = data->GetEntries();

  this->nnllblocks = 256;
  this->nllblocksize = 16;
  this->nnllthreads = this->nblocks * this->blocksize;

  this->expectations = new hemi::Array<float>(this->nsignals, true);
  for (size_t i=0; i<this->nsignals; i++) {
    this->expectations->writeOnlyHostPtr()[i] = signals[i].nexpected;
  }

  this->constraints = new hemi::Array<float>(this->nsignals, true);
  for (size_t i=0; i<this->nsignals; i++) {
    this->constraints->writeOnlyHostPtr()[i] = signals[i].constraint;
  }

  // Pj(xi) lookup table
  this->lut = MCMC::build_lut(signals, data);

  // if compiling device code, initialize the RNGs
#ifdef __CUDACC__
  this->rngs = new hemi::Array<RNGState>(this->nsignals, true);
  this->rngs->writeOnlyHostPtr();
  int bs = 256;
  int nb = this->nsignals / bs + 1;
  init_device_rngs<<<nb, bs>>>(this->nsignals, 1234, this->rngs->ptr());
#else
  this->rngs = static_cast<hemi::Array<RNGState>*>(NULL);
#endif
}


MCMC::~MCMC() {
  delete expectations;
  delete constraints;
  delete lut;
  delete rngs;
}


TNtuple* MCMC::operator()(unsigned nsteps, float burnin_fraction,
                    unsigned sync_interval) {
  int bs = 256;
  int nb = this->nsignals / bs + 1;

  float sigma = 0.125;  // FIXME make adaptive

  // buffers for current and proposed parameter vectors
  hemi::Array<float> current_vector(this->nsignals, true);
  for (size_t i=0; i<this->nsignals; i++) {
    current_vector.writeOnlyHostPtr()[i] = this->expectations->readOnlyHostPtr()[i];
  }

  hemi::Array<float> proposed_vector(this->nsignals, true);
  proposed_vector.writeOnlyHostPtr();  // touch to set valid

  // buffers for nll values at current and proposed parameter vectors
  hemi::Array<float> current_nll(1, true);
  current_nll.writeOnlyHostPtr();

  hemi::Array<float> proposed_nll(1, true);
  proposed_nll.writeOnlyHostPtr();

  // buffers for computing event term in nll
  hemi::Array<double> event_partial_sums(nnllthreads, true);
  hemi::Array<double> event_total_sum(1, true);    

  // calculate nll with initial parameters
  nll(current_vector, current_nll, event_partial_sums, event_total_sum);

  for (unsigned i=0; i<nsteps; i++) {
    HEMI_KERNEL_LAUNCH(pick_new_vector, nb, bs, 0, 0,
                       this->nsignals, this->rngs->ptr(), sigma,
                       current_vector.readOnlyPtr(),
                       proposed_vector.writeOnlyPtr());

    nll(proposed_vector, proposed_nll, event_partial_sums, event_total_sum);

    // compare, accept/reject, store

    std::cout << "current: " << current_nll.readOnlyHostPtr()[0] << ", "
              << "proposed: " << proposed_nll.readOnlyHostPtr()[0]
              << std::endl;

    if (i % sync_interval == 0 || i == nsteps - 1) {
      // dump gpu step buffer to TNtuple
    }
  }

  return static_cast<TNtuple*>(NULL);
}


void MCMC::nll(hemi::Array<float>& v, hemi::Array<float>& nll,
               hemi::Array<double>& event_partial_sums,
               hemi::Array<double>& event_total_sum) {
  // partial sums of event term
  HEMI_KERNEL_LAUNCH(nll_event_chunks, 16, 256, 0, 0,
                     this->lut->readOnlyPtr(),
                     v.readOnlyPtr(),
                     this->nevents, this->nsignals,
                     event_partial_sums.writeOnlyPtr());

  // total of event term
  HEMI_KERNEL_LAUNCH(nll_event_reduce, 1, this->nnllthreads, 0, 0,
                     this->nnllthreads,
                     event_partial_sums.readOnlyPtr(),
                     event_total_sum.writeOnlyPtr());

  // constraints + event term
  HEMI_KERNEL_LAUNCH(nll_total, 1, 1, 0, 0,
                     this->nsignals,
                     v.readOnlyPtr(),
                     this->expectations->readOnlyPtr(),
                     this->constraints->readOnlyPtr(),
                     event_total_sum.readOnlyPtr(),
                     nll.writeOnlyPtr());
}
    

hemi::Array<float>* MCMC::build_lut(const std::vector<Signal>& signals,
                                   TNtuple* data) {
  std::cout << "MCMC::build_lut: Building P(x) lookup table" << std::endl;
  int nevents = data->GetEntries();
  hemi::Array<float>* lut = new hemi::Array<float>(signals.size() * nevents,
                                                   true);

  std::vector<float> minima;
  for (size_t i=0; i<signals.size(); i++) {
    minima.push_back(signals[i].histogram->GetMinimum(0) * 0.0001);
  }

  float e;
  float r;
  data->SetBranchAddress("e", &e);
  data->SetBranchAddress("r", &r);

  for (int i=0; i<nevents; i++) {
    data->GetEntry(i);
    for (size_t j=0; j<signals.size(); j++) {
      double v = 0;
      if (signals[j].histogram->IsA() == TH2F::Class()) {
        v = dynamic_cast<TH2F*>(signals[j].histogram)->Interpolate(r, e);
      }
      else if (signals[j].histogram->IsA() == TH1D::Class()) {
        v = dynamic_cast<TH1D*>(signals[j].histogram)->Interpolate(e);
      }
      else {
        std::cerr << "build_lut: Unknown histogram class "
                  << signals[j].histogram->ClassName() << std::endl;
        assert(false);
      }

      if (v <= 0) {
        v = minima[j];
      }
      lut->writeOnlyHostPtr()[i * signals.size() + j] = v;
    }
  }

  return lut;
}

