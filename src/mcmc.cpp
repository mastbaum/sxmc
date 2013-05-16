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

  float* norms = new float[signals.size()];
  float* newnorms = new float[signals.size()];
  float* vals = new float[signals.size() + 1];  // norms + likelihood
  for (size_t i=0; i<signals.size(); i++) {
    norms[i] = signals[i].nexpected; // gRandom->Gaus(signals[i].nexpected, 5);
  }

  float sigma = 0.125; //0.75;  // FIXME: make this adaptive during burn-in, or settable

  NLL nll(signals, data);

  size_t naccepted = 0;
  size_t ntotal = 0;
  std::cout << "mcmc: Starting..." << std::endl;
  while (naccepted < nsteps) {
    if (ntotal % 1000 == 0) {
      std::cout << "mcmc: Accepted " << naccepted << " / " << nsteps
                << " (" << 100.0 * naccepted / ntotal << "\%, sigma = "
                << sigma << ")" << std::endl;
    }
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

      if (naccepted >= burnin_steps) {
        nt->Fill(vals);
      }
      //else {
      //  double accept_fraction = 1.0 * naccepted / ntotal;
      //  if (accept_fraction > 0.3) {
      //    sigma *= 1.01;
      //  }
      //  else if (accept_fraction < 0.25) {
      //    sigma *= 0.99;
      //  }
      //}

      naccepted++;
    }

    ntotal++;
  }

  delete[] norms;
  delete[] newnorms;
  delete[] vals;

  return nt;
}

