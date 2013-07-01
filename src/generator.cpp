#include <iostream>
#include <vector>
#include <assert.h>
#include <TRandom.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include "generator.h"
#include "signals.h"

std::vector<float> make_fake_dataset(std::vector<Signal> signals,
                                     std::vector<Systematic> systematics,
                                     std::vector<Observable> observables,
                                     float* params, bool poisson) {
  std::cout << "FakeDataGenerator::make_dataset: Generating dataset..."
            << std::endl;

  assert(observables.size() <= 3);

  std::vector<unsigned> observed(signals.size());

  size_t nevents = 0;
  for (size_t i=0; i<signals.size(); i++) {
    if (poisson) {
      observed[i] = gRandom->Poisson(params[i]);
    }
    else {
      observed[i] = TMath::Nint(params[i]);  // round to nearest integer
    }
    nevents += observed[i];
  }

  std::vector<float> data(nevents * observables.size());

  // evaluate the histogram once to force binning, and extract TH1 histogram
  std::vector<float> eval_points(signals.size());
  for (size_t i=0; i<signals.size(); i++) {
    eval_points[i] = 0;
  }

  hemi::Array<float> param_buffer(signals.size() + systematics.size(), true);
  for (size_t i=signals.size(); i<signals.size()+systematics.size(); i++) {
    param_buffer.writeOnlyHostPtr()[i] = params[i];
  }

  hemi::Array<float> value_buffer(1, true);
  hemi::Array<unsigned> norms_buffer(1, true);

  std::vector<TH1*> histograms(signals.size());

  for (size_t i=0; i<signals.size(); i++) {
    pdfz::Eval* p = signals[i].histogram;
    p->SetEvalPoints(eval_points);
    p->SetPDFValueBuffer(&value_buffer);
    p->SetNormalizationBuffer(&norms_buffer);
    p->SetParameterBuffer(&param_buffer);
    p->EvalAsync();
    p->EvalFinished();
    histograms[i] = dynamic_cast<pdfz::EvalHist*>(p)->CreateHistogram();
  }

  // generate event array by sampling ROOT histograms, including only events
  // that pass cuts
  std::vector<float> events(nevents * signals.size());
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

        events[i] = obs;
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

        events[i * signals.size() + 0] = obs0;
        events[i * signals.size() + 1] = obs1;
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

        events[i * signals.size() + 0] = obs0;
        events[i * signals.size() + 1] = obs1;
        events[i * signals.size() + 2] = obs2;
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

