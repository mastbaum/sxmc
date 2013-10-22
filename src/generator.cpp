#include <iostream>
#include <vector>
#include <assert.h>
#include <TRandom.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>

#include <sxmc/generator.h>
#include <sxmc/signals.h>

std::vector<float> make_fake_dataset(std::vector<Signal>& signals,
                                     std::vector<Systematic>& systematics,
                                     std::vector<Observable>& observables,
                                     std::vector<float> params, bool poisson) {
  std::cout << "make_fake_dataset: Generating dataset..." << std::endl;

  assert(observables.size() <= 3);

  // extract TH1 histogram from pdfz::Eval
  hemi::Array<double> param_buffer(signals.size() + systematics.size(), true);
  for (size_t i=0; i<signals.size() + systematics.size(); i++) {
    param_buffer.writeOnlyHostPtr()[i] = params[i];
  }

  hemi::Array<unsigned> norms_buffer(signals.size(), true);
  norms_buffer.writeOnlyHostPtr();

  std::vector<TH1*> histograms(signals.size());
  for (size_t i=0; i<signals.size(); i++) {
    pdfz::Eval* p = signals[i].histogram;
    p->SetNormalizationBuffer(&norms_buffer, i);
    p->SetParameterBuffer(&param_buffer, signals.size());
    histograms[i] = dynamic_cast<pdfz::EvalHist*>(p)->CreateHistogram();
  }

  std::vector<unsigned> observed(signals.size());
  size_t nevents = 0;
  for (size_t i=0; i<signals.size(); i++) {
    float r = params[i];

    if (poisson) {
      observed[i] = gRandom->Poisson(r);
    }
    else {
      observed[i] = TMath::Nint(r);  // round to the nearest integer
    }

    nevents += observed[i];
  }

  // generate event array by sampling ROOT histograms, including only events
  // that pass cuts
  size_t obs_id = 0;
  std::vector<float> events(nevents * observables.size());
  for (size_t i=0; i<signals.size(); i++) {
    std::cout << "make_fake_dataset: " << signals[i].name << ": "
              << observed[i] << " events (" << signals[i].nexpected
              << " expected)" << std::endl;

    if (histograms[i]->IsA() == TH1D::Class()) {
      TH1D* ht = dynamic_cast<TH1D*>(histograms[i]);

      double obs;
      for (unsigned j=0; j<observed[i]; j++) {
        do {
          obs = ht->GetRandom();
        } while(obs > observables[0].upper || obs < observables[0].lower);

        events[obs_id++] = obs;
      }
    }
    else if (histograms[i]->IsA() == TH2D::Class()) {
      TH2D* ht = dynamic_cast<TH2D*>(histograms[i]);

      double obs0;
      double obs1;
      for (unsigned j=0; j<observed[i]; j++) {
        do {
          ht->GetRandom2(obs0, obs1);
        } while(obs0 > observables[0].upper || obs0 < observables[0].lower ||
                obs1 > observables[1].upper || obs1 < observables[1].lower);

        events[obs_id++] = obs0;
        events[obs_id++] = obs1;
      }
    }
    else if (histograms[i]->IsA() == TH3D::Class()) {
      TH3D* ht = dynamic_cast<TH3D*>(histograms[i]);

      double obs0;
      double obs1;
      double obs2;
      for (unsigned j=0; j<observed[i]; j++) {
        do {
          ht->GetRandom3(obs0, obs1, obs2);
        } while(obs0 > observables[0].upper || obs0 < observables[0].lower ||
                obs1 > observables[1].upper || obs1 < observables[1].lower ||
                obs2 > observables[2].upper || obs2 < observables[2].lower);

        events[obs_id++] = obs0;
        events[obs_id++] = obs1;
        events[obs_id++] = obs2;
      }
    }
    else {
      std::cerr << "make_fake_dataset: Unknown histogram class: "
                << histograms[i]->ClassName() << std::endl;
      assert(false);
    }
  }

  return events;
}

