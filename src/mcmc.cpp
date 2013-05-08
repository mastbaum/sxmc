#include <iostream>
#include <vector>
#include <assert.h>
#include <TNtuple.h>
#include <TRandom.h>
#include <TMath.h>
#include "nll.h"
#include "signals.h"

TNtuple* mcmc(std::vector<Signal> signals, TNtuple* data, unsigned nsteps,
              float burnin_fraction) {
  unsigned burnin_steps = nsteps * burnin_fraction;

  // output ntuple
  std::string varlist;
  for (size_t i=0; i<signals.size(); i++) {
    varlist += (signals[i].name + ":");
  }
  varlist += "likelihood";

  TNtuple* nt = new TNtuple("likelihood", "Likelihood space", varlist.c_str());

  double* norms = new double[signals.size()];
  double* newnorms = new double[signals.size()];
  float* vals = new float[signals.size() + 1];  // norms + likelihood
  for (size_t i=0; i<signals.size(); i++) {
    norms[i] = gRandom->Gaus(signals[i].nexpected, 5);
  }

  float sigma = 5;  // FIXME: make this adaptive during burn-in, or settable

  NLL nll(signals, data);

  for (unsigned i=0; i<nsteps; i++) {
    double pcurrent = nll(norms);

    // jump
    for (size_t j=0; j<signals.size(); j++) {
      newnorms[j] = gRandom->Gaus(norms[j], sigma);
    }

    double pnew = nll(newnorms);

    // metropolis algorithm
    if (pnew < pcurrent || gRandom->Uniform() <= TMath::Exp(pcurrent - pnew)) {
      for (size_t j=0; j<signals.size(); j++) {
        assert(norms[j] == norms[j]);  // catch nans
        norms[j] = newnorms[j];
        vals[j] = norms[j];
      }
      vals[signals.size()] = pnew;
      if (i >= burnin_steps) {
        nt->Fill(vals);
      }
    }
  }

  delete[] norms;
  delete[] newnorms;
  delete[] vals;

  return nt;
}

