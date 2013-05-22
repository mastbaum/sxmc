#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <assert.h>
#include <hemi/hemi.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TNtuple.h>
#include <TRandom.h>
#include <TStopwatch.h>
#include <TDirectory.h>
#include "mcmc.h"
#include "signals.h"

#ifndef __HEMI_ARRAY_H__
#define __HEMI_ARRAY_H__
#include <hemi/array.h>
#endif

MCMC::MCMC(const std::vector<Signal>& signals, TNtuple* data) {
  this->nsignals = signals.size();
  this->nevents = data->GetEntries();

#ifdef __CUDACC__
  this->nnllblocks = 64;
  this->nllblocksize = 16;
#else
  this->nnllblocks = 1;
  this->nllblocksize = 1;
#endif
  this->nnllthreads = this->nnllblocks * this->nllblocksize;
  this->nreducethreads = 128;

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

  // list of variables for output ntuple
  for (size_t i=0; i<signals.size(); i++) {
    this->varlist += (signals[i].name + ":");
    this->signal_names.push_back(signals[i].name);
  }
  this->varlist += "likelihood";
  this->signal_names.push_back("likelihood");

  this->rngs = new hemi::Array<RNGState>(this->nsignals, true);

  // if compiling device code, initialize the RNGs
#ifdef __CUDACC__
  this->rngs->writeOnlyHostPtr();
  int bs = 128;
  int nb = this->nsignals / bs + 1;
  assert(nb < 8);
  init_device_rngs<<<nb, bs>>>(this->nsignals, 1234, this->rngs->ptr());
#else
  this->rngs->writeOnlyHostPtr();
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
  int bs = 128;
  int nb = this->nsignals / bs + 1;
  assert(nb < 8);

  unsigned burnin_steps = nsteps * burnin_fraction;

  // Ntuple to hold likelihood space
  TNtuple* nt = new TNtuple("lspace", "Likelihood space",
                            this->varlist.c_str());

  // buffers for current and proposed parameter vectors
  hemi::Array<float> current_vector(this->nsignals, true);
  for (size_t i=0; i<this->nsignals; i++) {
    current_vector.writeOnlyHostPtr()[i] = this->expectations->readOnlyHostPtr()[i];
  }

  hemi::Array<float> proposed_vector(this->nsignals, true);
  proposed_vector.writeOnlyHostPtr();  // touch to set valid

  // buffers for nll values at current and proposed parameter vectors
  hemi::Array<double> current_nll(1, true);
  current_nll.writeOnlyHostPtr();

  hemi::Array<double> proposed_nll(1, true);
  proposed_nll.writeOnlyHostPtr();

  // initial standard deviations for each dimension
  hemi::Array<float> sigma(this->nsignals, true);
  for (size_t i=0; i<this->nsignals; i++) {
    float expectation = this->expectations->readOnlyHostPtr()[i];
    float constraint = this->constraints->readOnlyHostPtr()[i];
    float width = (constraint > 0 ? constraint : 1.0 / sqrt(expectation));
    sigma.writeOnlyHostPtr()[i] = 5; //0.5 * width;
  }

  // buffers for computing event term in nll
  hemi::Array<double> event_partial_sums(this->nnllthreads, true);
  event_partial_sums.writeOnlyHostPtr();

  hemi::Array<double> event_total_sum(1, true);    
  event_total_sum.writeOnlyHostPtr();

  // buffer of accepted jumps, transferred from gpu periodically
  hemi::Array<int> jump_counter(1, true);
  jump_counter.writeOnlyHostPtr()[0] = 0;

  hemi::Array<float> jump_buffer(sync_interval * (this->nsignals + 1), true);
  float* jump_vector = new float[this->nsignals + 1];

  // calculate nll with initial parameters
  const float* cv = current_vector.readOnlyHostPtr();
  nll(current_vector.readOnlyPtr(), current_nll.writeOnlyPtr(),
      event_partial_sums.ptr(), event_total_sum.ptr());

  // perform random walk
  TStopwatch timer;
  timer.Start();
  for (unsigned i=0; i<nsteps; i++) {
    // re-tune jump distribution based on burn-in phase
    if (i == burnin_steps) {
      std::cout << "MCMC: Burn-in phase completed after " << burnin_steps
                << " steps" << std::endl;

      // fit a Gaussian in each dimension to estimate distribution width
      for (size_t j=0; j<this->nsignals; j++) {
        std::string name = this->signal_names[j];
        nt->Draw((name + ">>hsproj").c_str());
        TH1F* hsproj = (TH1F*) gDirectory->Get("hsproj");

        if (!hsproj) {
          std::cerr << "MCMC: failed to get signal projection" << std::endl;
          continue;
        }

        hsproj->Fit("gaus", "q");
        TF1* fsproj = hsproj->GetFunction("gaus");

        if (!fsproj) {
          std::cerr << "MCMC: failed to fit signal projection" << std::endl;
          continue;
        }

        double fit_width = fsproj->GetParameter(2);
        const float scale_factor = 2.6;  // e!?

        sigma.writeOnlyHostPtr()[j] = scale_factor * fit_width;
        std::cout << "MCMC: Rescaling jump sigma: " << name << " "
                  << scale_factor * fit_width << std::endl;

        hsproj->Delete();
      }

      nt->Reset();
    }

    HEMI_KERNEL_LAUNCH(pick_new_vector, 8, 8, 0, 0,
                       this->nsignals, this->rngs->ptr(),
                       sigma.readOnlyPtr(),
                       current_vector.readOnlyPtr(),
                       proposed_vector.writeOnlyPtr());

    nll(proposed_vector.readOnlyPtr(), proposed_nll.writeOnlyPtr(),
        event_partial_sums.ptr(), event_total_sum.ptr());

    //std::cout << "current: " << current_nll.readOnlyHostPtr()[0] << ", "
    //          << "proposed: " << proposed_nll.readOnlyHostPtr()[0]
    //          << std::endl;

    HEMI_KERNEL_LAUNCH(jump_decider, 1, 1, 0, 0,
                       this->rngs->ptr(),
                       current_nll.ptr(),
                       proposed_nll.readOnlyPtr(),
                       current_vector.ptr(),
                       proposed_vector.readOnlyPtr(),
                       this->nsignals,
                       jump_counter.ptr(),
                       jump_buffer.writeOnlyPtr());

    // flush the jump buffer periodically
    if (i % sync_interval == 0 || i == nsteps - 1) {
      int njumps = jump_counter.readOnlyHostPtr()[0];
      std::cout << i << " " << sync_interval << " " << nsteps << ": JUMPS IN BUFFER: " << njumps << std::endl;
      for (int j=0; j<njumps; j++) {
         // first nsignals elements are normalizations
         for (unsigned k=0; k<this->nsignals; k++) {
           int idx = j * (this->nsignals + 1) + k;
           jump_vector[k] = jump_buffer.readOnlyHostPtr()[idx];
         }
         // last element is the likelihood
         jump_vector[this->nsignals] = jump_buffer.readOnlyHostPtr()[j*(this->nsignals+1)+this->nsignals];

         nt->Fill(jump_vector);
      }

      // reset counter
      jump_counter.writeOnlyHostPtr()[0] = 0;
    }
  }

  std::cout << "MCMC: Elapsed time: " << timer.RealTime() << std::endl;

  delete[] jump_vector;

  return nt;
}


void MCMC::nll(const float* v, double* nll,
               double* event_partial_sums,
               double* event_total_sum) {
  // partial sums of event term
  HEMI_KERNEL_LAUNCH(nll_event_chunks, this->nnllblocks,
                     this->nllblocksize, 0, 0,
                     this->lut->readOnlyPtr(),
                     v,
                     this->nevents, this->nsignals,
                     event_partial_sums);

  // total of event term
  HEMI_KERNEL_LAUNCH(nll_event_reduce, 1, this->nreducethreads,  
		     this->nreducethreads * sizeof(double) /* shared memory */, 0,
                     this->nnllthreads,
                     event_partial_sums,
                     event_total_sum);

  // constraints + event term
  HEMI_KERNEL_LAUNCH(nll_total, 1, 1, 0, 0,
                     this->nsignals,
                     v,
                     this->expectations->readOnlyPtr(),
                     this->constraints->readOnlyPtr(),
                     event_total_sum,
                     nll);
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

