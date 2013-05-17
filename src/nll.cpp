#include <iostream>
#include <vector>
#include <cmath>
#include <assert.h>
#include <cuda.h>
#include <hemi/hemi.h>
#ifndef __HEMI_ARRAY_H__
#define __HEMI_ARRAY_H__
#include <hemi/array.h>
#endif
#include <TNtuple.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH1.h>
#include "signals.h"
#include "nll.h"

#define DEBUG  // let HEMI print errors

HEMI_KERNEL(ll)(const float* lut, const float* pars, const size_t ne,
                const size_t ns, double* sums) {
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


NLL::NLL(const std::vector<Signal>& signals, TNtuple* data) {
  this->nsignals = signals.size();
  this->nevents = data->GetEntries();

  this->expectations = new hemi::Array<float>(nsignals, true);
  this->constraints = new hemi::Array<float>(nsignals, true);

  for (size_t i=0; i<nsignals; i++) {
    this->expectations->writeOnlyHostPtr()[i] = signals[i].nexpected;
    this->constraints->writeOnlyHostPtr()[i] = signals[i].constraint;
  }

  this->lut = NLL::build_lut(signals, data);

  // pre-allocate buffers for the normalizations and output sums,
  // which change on every call
  this->normalizations = new hemi::Array<float>(nsignals, true);

  this->blocksize = 256;
  this->nblocks = 16;
  this->nthreads = this->nblocks * this->blocksize;
  this->sums = new hemi::Array<double>(this->nthreads, true);
  this->sums->writeOnlyHostPtr();  // no-op to set hemi::Array validity
}


NLL::~NLL() {
  delete this->lut;
  delete this->expectations;
  delete this->constraints;
  delete this->normalizations;
  delete this->sums;
}


hemi::Array<float>* NLL::build_lut(const std::vector<Signal>& signals, TNtuple* data) {
  std::cout << "NLL::build_lut: Building P(x) lookup table" << std::endl;
  int nevents = data->GetEntries();
  hemi::Array<float>* lut = new hemi::Array<float>(signals.size() * nevents, true);

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


double NLL::operator()(float* norms) {
  double result = 0;

  // N + fractional constraints
  for (size_t i=0; i<this->nsignals; i++) {
    if (norms[i] < 0) {
      return 1e10;
    }
    result += norms[i];
    if (this->constraints->readOnlyHostPtr()[i] > 0) {
      result += 0.5 * pow((norms[i]/this->expectations->readOnlyHostPtr()[i] - 1) /
                          this->constraints->readOnlyHostPtr()[i], 2);
    }
  }

  for (size_t i=0; i<this->nsignals; i++) {
    this->normalizations->writeOnlyHostPtr()[i] = norms[i];
  }

  HEMI_KERNEL_LAUNCH(ll, this->nblocks, this->blocksize, 0, 0,
                     this->lut->readOnlyPtr(),
                     this->normalizations->readOnlyPtr(),
                     this->nevents, this->nsignals,
                     this->sums->ptr());

  double sum = 0;
  for (size_t i=0; i<this->nthreads; i++) {
    sum += this->sums->readOnlyHostPtr()[i];
  }

  result -= sum;

  return result;
}

