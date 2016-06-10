#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <cassert>
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
#include <sxmc/signal.h>
#include <sxmc/observable.h>
#include <sxmc/systematic.h>
#include <sxmc/likelihood.h>

#ifndef __HEMI_ARRAY_H__
#define __HEMI_ARRAY_H__
#include <hemi/array.h>
#endif

MCMC::MCMC(const std::vector<Source>& sources,
           const std::vector<Signal>& signals,
           const std::vector<Systematic>& systematics,
           const std::vector<Observable>& observables) {
  this->nsources = sources.size();
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

  // Total number of systematic parameters
  size_t npars = 0;
  for (size_t i=0; i<systematics.size(); i++) {
    npars += systematics[i].npars;
  }

  // Set mean/expectation and sigma for all fit parameters
  this->nparameters = this->nsources + npars;
  this->parameter_means = new hemi::Array<double>(this->nparameters, true);
  this->parameter_sigma = new hemi::Array<double>(this->nparameters, true);
  this->parameter_fixed.resize(this->nparameters);
  this->nfloat = 0; 

  for (size_t i=0; i<this->nsources; i++) {
    this->parameter_means->writeOnlyHostPtr()[i] = sources[i].mean;
    this->parameter_sigma->writeOnlyHostPtr()[i] = sources[i].sigma;
    this->parameter_fixed[i] = sources[i].fixed;
    this->nfloat += (sources[i].fixed ? 0 : 1);
  }

  this->systematics_fixed = true;

  size_t k = this->nsources;
  for (size_t i=0; i<systematics.size(); i++) {
    if (!systematics[i].fixed) {
      this->systematics_fixed = false;
    }
    for (size_t j=0; j<systematics[i].npars; j++) {
      this->parameter_means->writeOnlyHostPtr()[k] = systematics[i].means[j];
      this->parameter_sigma->writeOnlyHostPtr()[k] = systematics[i].sigmas[j];
      this->parameter_fixed[k] = systematics[i].fixed;
      this->nfloat += (systematics[i].fixed ? 0 : 1);
      k++;
    }
  }

  if (systematics.size() > 0 && this->systematics_fixed) {
    std::cout << "MCMC::MCMC: All systematics are fixed, will not reevaluate "
              << "PDFs." << std::endl;
  }

  // Signal expectations, datasets, pdfz::Eval histograms
  this->pdfs.resize(this->nsignals);
  this->nexpected = new hemi::Array<double>(this->nsignals, true);
  this->n_mc = new hemi::Array<unsigned>(this->nsignals, true);
  this->source_id = new hemi::Array<short>(this->nsignals, true);
  for (size_t i=0; i<this->nsignals; i++) {
    this->pdfs[i] = signals[i].histogram;
    this->nexpected->writeOnlyHostPtr()[i] = signals[i].nexpected;
    this->n_mc->writeOnlyHostPtr()[i] = signals[i].n_mc;
    this->source_id->writeOnlyHostPtr()[i] = signals[i].source.index;
  }

  // List of parameters for output ntuple
  for (size_t i=0; i<sources.size(); i++) {
    this->varlist += (sources[i].name + ":");
    this->parameter_names.push_back(sources[i].name);
  }
  for (size_t i=0; i<systematics.size(); i++) {
    for (size_t j=0; j<systematics[i].npars; j++) {
      std::ostringstream oss;
      oss << systematics[i].name << "_" << j;
      this->varlist += (oss.str() + ":");
      this->parameter_names.push_back(oss.str());
    }
  }
  this->varlist += "likelihood";
  this->parameter_names.push_back("likelihood");

  this->rngs = new hemi::Array<RNGState>(this->nparameters, true);

  // If compiling device code, initialize the RNGs
#ifdef __CUDACC__
  this->rngs->writeOnlyHostPtr();
  int bs = 128;
  int nb = this->nparameters / bs + 1;
  assert(nb < 8);
  init_device_rngs<<<nb, bs>>>(this->nparameters,
                               gRandom->GetSeed(),  // Random seed
                               this->rngs->ptr());
#else
  this->rngs->writeOnlyHostPtr();
#endif
}


MCMC::~MCMC() {
  delete parameter_means;
  delete parameter_sigma;
  delete nexpected;
  delete n_mc;
  delete source_id;
  delete rngs;
}


LikelihoodSpace*
MCMC::operator()(std::vector<float>& data, unsigned nsteps,
                 float burnin_fraction, const bool debug_mode,
                 unsigned sync_interval) {
  // CUDA/hemi block sizes
  int bs = 128;
  int nb = this->nsignals / bs + 1;
  assert(nb < 8);

  unsigned burnin_steps = nsteps * burnin_fraction;

  // Ntuple to hold likelihood space
  TNtuple* nt = new TNtuple("lspace", "Likelihood space",
                            this->varlist.c_str());

  // Buffers for current and proposed parameter vectors
  hemi::Array<double> current_vector(this->nparameters, true);
  for (size_t i=0; i<this->nparameters; i++) {
    current_vector.writeOnlyHostPtr()[i] = \
      this->parameter_means->readOnlyHostPtr()[i];
  }

  hemi::Array<double> proposed_vector(this->nparameters, true);
  proposed_vector.writeOnlyHostPtr();  // Touch to set valid

  // Buffer for normalizations after application of systematics
  hemi::Array<unsigned> normalizations(this->nsignals, true);
  normalizations.writeOnlyHostPtr();

  // Buffers for computing event term in nll
  hemi::Array<double> event_partial_sums(this->nnllthreads, true);
  event_partial_sums.writeOnlyHostPtr();

  hemi::Array<double> event_total_sum(1, true);
  event_total_sum.writeOnlyHostPtr();

  // Buffer of jumps, transferred from device memory periodically
  hemi::Array<int> jump_counter(1, true);
  jump_counter.writeOnlyHostPtr()[0] = 0;

  hemi::Array<int> accept_counter(1, true);
  accept_counter.writeOnlyHostPtr()[0] = 0;

  hemi::Array<float> jump_buffer(sync_interval * (this->nparameters+1), true);

  float* jump_vector = new float[this->nparameters + 1];

  // Buffers for nll values at current and proposed parameter vectors
  hemi::Array<double> current_nll(1, true);
  current_nll.writeOnlyHostPtr();

  hemi::Array<double> proposed_nll(1, true);
  proposed_nll.writeOnlyHostPtr();

  // Set initial proposal distribution widths for each dimension
  hemi::Array<float> jump_width(this->nparameters, true);
  const float scale_factor = 2.4 * 2.4 / this->nfloat;  // Haario, 2001
  std::cout << "MCMC: Intial jump sigma" << std::endl;
  for (size_t i=0; i<this->nparameters; i++) {
    std::cout << " " << this->parameter_names[i] << ": ";

    if (this->parameter_fixed[i]) {
      std::cout << "fixed" << std::endl;
      jump_width.writeOnlyHostPtr()[i] = -1;
      continue;
    }

    float mean = this->parameter_means->readOnlyHostPtr()[i];
    float sigma = this->parameter_sigma->readOnlyHostPtr()[i];
    float width = 0.1;

    if (sigma > 0) {
      width = sigma;
    }
    else if (i < nsignals) {
      float m = std::max(mean, (float) 10);
      width = sqrt(m) / m;
    }
    else {
      width = sqrt(std::max(mean, (float) 1));
    }

    jump_width.writeOnlyHostPtr()[i] = 0.1 * width * scale_factor;

    std::cout << jump_width.readOnlyHostPtr()[i] << std::endl;
  }

  // Set up histogram and perform initial evaluation
  size_t nevents = data.size() / (this->nobservables + 1);
  hemi::Array<float> lut(nevents * this->nsignals, true);
  for (size_t i=0; i<this->pdfs.size(); i++) {
    pdfz::Eval* p = this->pdfs[i];
    p->SetEvalPoints(data);
    p->SetPDFValueBuffer(&lut, i * nevents, 1);
    p->SetNormalizationBuffer(&normalizations, i);
    p->SetParameterBuffer(&current_vector, this->nsources);
    p->EvalAsync();
    p->EvalFinished();
    p->SetParameterBuffer(&proposed_vector, this->nsources);
  }

  // Calculate NLL with initial parameters
  nll(lut.readOnlyPtr(), nevents,
      current_vector.readOnlyPtr(), current_nll.writeOnlyPtr(),
      this->nexpected->readOnlyPtr(), this->n_mc->readOnlyPtr(),
      this->source_id->readOnlyPtr(),
      normalizations.readOnlyPtr(),
      event_partial_sums.ptr(), event_total_sum.ptr());

  HEMI_KERNEL_LAUNCH(pick_new_vector, 1, 64, 0, 0,
                     this->nparameters, this->rngs->ptr(),
                     jump_width.readOnlyPtr(),
                     current_vector.readOnlyPtr(),
                     proposed_vector.writeOnlyPtr());

  // Perform random walk
  TStopwatch timer;
  timer.Start();
  for (unsigned i=0; i<nsteps; i++) {
    // If systematics are varying, re-evaluate the pdfs
    // At some point, look for ways of doing this less often?
    if (this->nsystematics > 0 && !this->systematics_fixed) {
      for (size_t j=0; j<this->pdfs.size(); j++) {
        pdfs[j]->EvalAsync();
      }
      for (size_t j=0; j<this->pdfs.size(); j++) {
        pdfs[j]->EvalFinished();
      }
    }

    // Re-tune jump distribution based on burn-in phase steps
    if (i == burnin_steps || i == 2 * burnin_steps) {
      std::cout << "MCMC: Burn-in phase completed after " << burnin_steps
                << " steps" << std::endl;

      // Rescale jumps in each dimension based on RMS during burn-in
      std::cout << "MCMC: Rescaling jump sigma" << std::endl;
      for (size_t j=0; j<this->nparameters; j++) {
        std::string name = this->parameter_names[j];

        if (this->parameter_fixed[j]) {
          std::cout << " " << name << ": fixed" << std::endl;
          continue;
        }

        nt->Draw((name + ">>hsproj").c_str());
        TH1F* hsproj = (TH1F*) gDirectory->Get("hsproj");
        assert(hsproj);

        double fit_width = (hsproj->GetRMS() > 0 ? hsproj->GetRMS() :
                            jump_width.readOnlyHostPtr()[j]);

        std::cout << " " << name << ": "
                  << jump_width.readOnlyHostPtr()[j] << " -> ";

        jump_width.writeOnlyHostPtr()[j] = scale_factor * fit_width;

        std::cout << jump_width.readOnlyHostPtr()[j] << std::endl;

        hsproj->Delete();
      }

      jump_width.writeOnlyHostPtr();

      // Save all steps when in debug mode
      if (!debug_mode) {
        nt->Reset();
      }
    }

    // Partial sums of event term
    HEMI_KERNEL_LAUNCH(nll_event_chunks, this->nnllblocks,
                       this->nllblocksize, 0, 0,
                       lut.readOnlyPtr(),
                       proposed_vector.readOnlyPtr(),
                       nevents, this->nsignals,
                       this->nexpected->readOnlyPtr(),
                       this->n_mc->readOnlyPtr(),
                       this->source_id->readOnlyPtr(),
                       normalizations.readOnlyPtr(),
                       event_partial_sums.ptr());

    // Accept/reject the jump, add current position to the buffer
    HEMI_KERNEL_LAUNCH(finish_nll_jump_pick_combo, 1, this->nreducethreads,
                       this->nreducethreads * sizeof(double), 0,
                       this->nnllthreads,
                       event_partial_sums.ptr(),
                       this->nsignals,
                       this->nsources,
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
                       this->nexpected->readOnlyPtr(),
                       this->n_mc->readOnlyPtr(),
                       this->source_id->readOnlyPtr(),
                       normalizations.readOnlyPtr(),
                       debug_mode);

    // Flush the jump buffer periodically
    if (i % sync_interval == 0 || i == nsteps - 1 ||
        i == burnin_steps - 1 || i == 2 * burnin_steps - 1) {
      int njumps = jump_counter.readOnlyHostPtr()[0];
      int naccepted = accept_counter.readOnlyHostPtr()[0];

      std::cout << "MCMC: Step " << i + 1 << "/" << nsteps
                << " (" << njumps << " in buffer, "
                << naccepted << " accepted)" << std::endl;

      for (int j=0; j<njumps; j++) {
         // First nparameters elements are parameters...
         for (unsigned k=0; k<this->nparameters; k++) {
           int idx = j * (this->nparameters + 1) + k;
           jump_vector[k] = jump_buffer.readOnlyHostPtr()[idx];
         }
         // ... Last element is the (negative log) likelihood
         jump_vector[this->nparameters] = \
           jump_buffer.readOnlyHostPtr()[j * (this->nparameters + 1) +
                                         this->nparameters];

         nt->Fill(jump_vector);
      }

      // Reset counters
      jump_counter.writeOnlyHostPtr()[0] = 0;
      accept_counter.writeOnlyHostPtr()[0] = 0;
    }
  }

  std::cout << "MCMC: Elapsed time: " << timer.RealTime() << std::endl;

  LikelihoodSpace* lspace = new LikelihoodSpace(nt);

  delete[] jump_vector;

  return lspace;
}


void MCMC::nll(const float* lut, size_t nevents, const double* v, double* nll,
               const double* nexpected, const unsigned* n_mc,
               const short* source_id,
               const unsigned* norms,
               double* event_partial_sums, double* event_total_sum) {
  // Partial sums of event term
  HEMI_KERNEL_LAUNCH(nll_event_chunks,
                     this->nnllblocks, this->nllblocksize, 0, 0,
                     lut, v, nevents, this->nsignals,
                     nexpected, n_mc, source_id,
                     norms,
                     event_partial_sums);

  // Total of event term
  HEMI_KERNEL_LAUNCH(nll_event_reduce, 1, this->nreducethreads,  
                     this->nreducethreads * sizeof(double), 0,
                     this->nnllthreads, event_partial_sums, event_total_sum);

  // Constraints + event term
  HEMI_KERNEL_LAUNCH(nll_total, 1, 1, 0, 0,
                     this->nparameters, v, this->nsignals, this->nsources,
                     this->parameter_means->readOnlyPtr(),
                     this->parameter_sigma->readOnlyPtr(),
                     event_total_sum, nexpected, n_mc, source_id,
                     norms, nll);
}

