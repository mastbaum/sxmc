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

#include <sxmc/mcmc.h>
#include <sxmc/signals.h>
#include <sxmc/likelihood.h>

#ifndef __HEMI_ARRAY_H__
#define __HEMI_ARRAY_H__
#include <hemi/array.h>
#endif

MCMC::MCMC(const std::vector<Signal>& signals,
           const std::vector<Systematic>& systematics,
           const std::vector<Observable>& observables) {
  this->nsignals = signals.size();
  this->nsystematics = systematics.size();
  this->nobservables = observables.size();

#ifdef __CUDACC__
  this->nnllblocks = 64;
  this->nllblocksize = 256;
#else
  this->nnllblocks = 1;
  this->nllblocksize = 1;
#endif
  this->nnllthreads = this->nnllblocks * this->nllblocksize;
  this->nreducethreads = 128;

  // set mean/expectation and sigma for all parameters
  this->nparameters = this->nsignals + this->nsystematics;
  this->parameter_means = new hemi::Array<double>(this->nparameters, true);
  this->parameter_sigma = new hemi::Array<double>(this->nparameters, true);
  for (size_t i=0; i<this->nsignals; i++) {
    this->parameter_means->writeOnlyHostPtr()[i] = signals[i].nexpected;
    this->parameter_sigma->writeOnlyHostPtr()[i] = signals[i].sigma;
  }
  for (size_t i=0; i<this->nsystematics; i++) {
    this->parameter_means->writeOnlyHostPtr()[this->nsignals + i] = \
      systematics[i].mean;
    this->parameter_sigma->writeOnlyHostPtr()[this->nsignals + i] = \
      systematics[i].sigma;
  }

  // references to pdfz::Eval histograms
  this->pdfs.resize(this->nsignals);
  for (size_t i=0; i<this->nsignals; i++) {
    this->pdfs[i] = signals[i].histogram;
  }

  // list of parameters for output ntuple
  for (size_t i=0; i<signals.size(); i++) {
    this->varlist += (signals[i].name + ":");
    this->parameter_names.push_back(signals[i].name);
  }
  for (size_t i=0; i<systematics.size(); i++) {
    this->varlist += (systematics[i].name + ":");
    this->parameter_names.push_back(systematics[i].name);
  }
  this->varlist += "likelihood";
  this->parameter_names.push_back("likelihood");

  this->rngs = new hemi::Array<RNGState>(this->nparameters, true);

  // if compiling device code, initialize the RNGs
#ifdef __CUDACC__
  this->rngs->writeOnlyHostPtr();
  int bs = 128;
  int nb = this->nparameters / bs + 1;
  assert(nb < 8);
  init_device_rngs<<<nb, bs>>>(this->nparameters, 1234, this->rngs->ptr());
#else
  this->rngs->writeOnlyHostPtr();
#endif
}


MCMC::~MCMC() {
  delete parameter_means;
  delete parameter_sigma;
  delete rngs;
}


LikelihoodSpace* MCMC::operator()(std::vector<float>& data, std::vector<int>& weights, unsigned nsteps,
                                  float burnin_fraction, const bool debug_mode,
                                  unsigned sync_interval) {
  // cuda/hemi block sizes
  int bs = 128;
  int nb = this->nsignals / bs + 1;
  assert(nb < 8);

  unsigned burnin_steps = nsteps * burnin_fraction;

  // Ntuple to hold likelihood space
  TNtuple* nt = new TNtuple("lspace", "Likelihood space",
                            this->varlist.c_str());

  // buffers for current and proposed parameter vectors
  hemi::Array<double> current_vector(this->nparameters, true);
  for (size_t i=0; i<this->nparameters; i++) {
    current_vector.writeOnlyHostPtr()[i] = \
      this->parameter_means->readOnlyHostPtr()[i];
  }

  hemi::Array<double> proposed_vector(this->nparameters, true);
  proposed_vector.writeOnlyHostPtr();  // touch to set valid

  // buffer for normalizations after application of systematics
  hemi::Array<unsigned> normalizations(this->nsignals, true);
  normalizations.writeOnlyHostPtr();

  // buffers for nll values at current and proposed parameter vectors
  hemi::Array<double> current_nll(1, true);
  current_nll.writeOnlyHostPtr();

  hemi::Array<double> proposed_nll(1, true);
  proposed_nll.writeOnlyHostPtr();

  // create hemi buffer for weighting data points
  hemi::Array<int> dataweights(data.size(), true);
  dataweights.copyFromHost(&weights.front(), weights.size());

  // initial standard deviations for each dimension
  hemi::Array<float> jump_width(this->nparameters, true);
  const float scale_factor = 2.4 * 2.4 / this->nparameters;  // Haario, 2001
  for (size_t i=0; i<this->nparameters; i++) {
    float mean = max(this->parameter_means->readOnlyHostPtr()[i], 10.0);
    float sigma = this->parameter_sigma->readOnlyHostPtr()[i];
    float width = (sigma > 0 ? sigma : sqrt(mean));
    jump_width.writeOnlyHostPtr()[i] = 0.1 * width * scale_factor;
  }

  // buffers for computing event term in nll
  hemi::Array<double> event_partial_sums(this->nnllthreads, true);
  event_partial_sums.writeOnlyHostPtr();

  hemi::Array<double> event_total_sum(1, true);    
  event_total_sum.writeOnlyHostPtr();

  // buffer of jumps, transferred from gpu periodically
  hemi::Array<int> jump_counter(1, true);
  jump_counter.writeOnlyHostPtr()[0] = 0;

  hemi::Array<int> accept_counter(1, true);
  accept_counter.writeOnlyHostPtr()[0] = 0;

  hemi::Array<float> jump_buffer(sync_interval * (this->nparameters + 1),
                                 true);

  float* jump_vector = new float[this->nparameters + 1];

  // set up histogram and perform initial evaluation
  size_t nevents = data.size() / this->nobservables;
  hemi::Array<float> lut(nevents * this->nsignals, true);
  for (size_t i=0; i<this->pdfs.size(); i++) {
    pdfz::Eval* p = this->pdfs[i];
    p->SetEvalPoints(data);
    p->SetPDFValueBuffer(&lut, i * nevents, 1);
    p->SetNormalizationBuffer(&normalizations, i);
    p->SetParameterBuffer(&current_vector, this->nsignals);
    p->EvalAsync();
    p->EvalFinished();
    p->SetParameterBuffer(&proposed_vector, this->nsignals);
  }

  // calculate nll with initial parameters
  nll(lut.readOnlyPtr(), dataweights.readOnlyPtr(), nevents,
      current_vector.readOnlyPtr(), current_nll.writeOnlyPtr(),
      event_partial_sums.ptr(), event_total_sum.ptr());

  HEMI_KERNEL_LAUNCH(pick_new_vector, 1, 64, 0, 0,
                     this->nparameters, this->rngs->ptr(),
                     jump_width.readOnlyPtr(),
                     current_vector.readOnlyPtr(),
                     proposed_vector.writeOnlyPtr());

  // perform random walk
  TStopwatch timer;
  timer.Start();
  for (unsigned i=0; i<nsteps; i++) {
    // if systematics are varying, re-evaluate the pdfs
    if (this->nsystematics > 0) {
      for (size_t i=0; i<this->pdfs.size(); i++) {
        pdfs[i]->EvalAsync();
      }
      for (size_t i=0; i<this->pdfs.size(); i++) {
        pdfs[i]->EvalFinished();
      }
    }

    // re-tune jump distribution based on burn-in phase
    if (i == burnin_steps || i == 2 * burnin_steps) {
      std::cout << "MCMC: Burn-in phase completed after " << burnin_steps
                << " steps" << std::endl;

      // rescale jumps in each dimension based on RMS during burn-in
      for (size_t j=0; j<this->nparameters; j++) {
        std::string name = this->parameter_names[j];
        nt->Draw((name + ">>hsproj").c_str());
        TH1F* hsproj = (TH1F*) gDirectory->Get("hsproj");

        double fit_width = hsproj->GetRMS();

        std::cout << "MCMC: Rescaling jump sigma: " << name << ": "
                  << jump_width.readOnlyHostPtr()[j] << " -> ";

        jump_width.writeOnlyHostPtr()[j] = scale_factor * fit_width;

        std::cout << jump_width.readOnlyHostPtr()[j] << std::endl;

        hsproj->Delete();
      }
      // save all steps when in debug mode
      if (!debug_mode) {
        nt->Reset();
      }
    }

    // partial sums of event term
    HEMI_KERNEL_LAUNCH(nll_event_chunks, this->nnllblocks,
                       this->nllblocksize, 0, 0,
                       lut.readOnlyPtr(), dataweights.readOnlyPtr(),
                       proposed_vector.readOnlyPtr(),
                       nevents, this->nsignals,
                       event_partial_sums.ptr());

    // accept/reject the jump, add current position to the buffer
    HEMI_KERNEL_LAUNCH(finish_nll_jump_pick_combo, 1, this->nreducethreads,
                       this->nreducethreads * sizeof(double), 0,
                       this->nnllthreads,
                       event_partial_sums.ptr(),
                       this->nsignals, 
                       this->parameter_means->readOnlyPtr(),
                       this->parameter_sigma->readOnlyPtr(),
                       this->rngs->ptr(),
                       current_nll.ptr(),
                       proposed_nll.ptr(),
                       current_vector.ptr(),
                       proposed_vector.ptr(),
                       accept_counter.ptr(),
                       jump_counter.ptr(),
                       jump_buffer.writeOnlyPtr(),
                       this->nparameters,
                       jump_width.readOnlyPtr(),
                       debug_mode);

    // flush the jump buffer periodically
    if (i % sync_interval == 0 || i == nsteps - 1 || i == burnin_steps - 1) {
      int njumps = jump_counter.readOnlyHostPtr()[0];
      int naccepted = accept_counter.readOnlyHostPtr()[0];
      std::cout << "MCMC: Step " << i << "/" << nsteps
                << " (" << njumps << " in buffer, "
                << naccepted << " accepted)" << std::endl;
      for (int j=0; j<njumps; j++) {
         // first nsignals elements are normalizations
         for (unsigned k=0; k<this->nparameters; k++) {
           int idx = j * (this->nparameters + 1) + k;
           jump_vector[k] = jump_buffer.readOnlyHostPtr()[idx];
         }
         // last element is the likelihood
         jump_vector[this->nparameters] = \
           jump_buffer.readOnlyHostPtr()[j * (this->nparameters + 1) +
                                         this->nparameters];

         nt->Fill(jump_vector);
      }

      // reset counters
      jump_counter.writeOnlyHostPtr()[0] = 0;
      accept_counter.writeOnlyHostPtr()[0] = 0;
    }
  }

  std::cout << "MCMC: Elapsed time: " << timer.RealTime() << std::endl;


  LikelihoodSpace* lspace = new LikelihoodSpace(nt);

  delete[] jump_vector;

  return lspace;
}


void MCMC::nll(const float* lut, const int* dataweights, size_t nevents, const double* v, double* nll,
               double* event_partial_sums, double* event_total_sum) {
  // partial sums of event term
  HEMI_KERNEL_LAUNCH(nll_event_chunks,
                     this->nnllblocks, this->nllblocksize, 0, 0,
                     lut, dataweights, v, nevents, this->nsignals, event_partial_sums);

  // total of event term
  HEMI_KERNEL_LAUNCH(nll_event_reduce, 1, this->nreducethreads,  
                     this->nreducethreads * sizeof(double), 0,
                     this->nnllthreads, event_partial_sums, event_total_sum);

  // constraints + event term
  HEMI_KERNEL_LAUNCH(nll_total, 1, 1, 0, 0,
                     this->nparameters, v, this->nsignals,
                     this->parameter_means->readOnlyPtr(),
                     this->parameter_sigma->readOnlyPtr(),
                     event_total_sum, nll);
}

