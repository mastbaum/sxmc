#include <iostream>
#include <utility>
#include <vector>
#include <assert.h>
#include <TRandom.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>

#include <sxmc/generator.h>
#include <sxmc/signals.h>


std::pair<std::vector<float>, std::vector<int> >
make_fake_dataset(std::vector<Signal>& signals,
                  std::vector<Systematic>& systematics,
                  std::vector<Observable>& observables,
                  std::vector<double> params, bool poisson, int maxsamples) {
  std::cout << "make_fake_dataset: Generating dataset..." << std::endl;

  std::vector<double> syst_vals;
  for (size_t i=signals.size(); i<signals.size() + systematics.size(); i++) {
    syst_vals.push_back(params[i]);
  }

  std::vector<float> upper;
  std::vector<float> lower;
  for (size_t i=0; i<observables.size(); i++) {
    upper.push_back(observables[i].upper);
    lower.push_back(observables[i].lower);
  }

  std::vector<float> events;
  std::vector<int> weights;
  std::vector<unsigned> observed(signals.size());
  for (size_t i=0; i<signals.size(); i++) {
    observed[i] = \
      dynamic_cast<pdfz::EvalHist*>(signals[i].histogram)->RandomSample(
        events, weights, params[i], syst_vals, upper, lower,
        poisson, maxsamples);

    std::cout << "make_fake_dataset: " << signals[i].name << ": "
              << observed[i] << " events (" << signals[i].nexpected
              << " expected)" << std::endl;
  }

  return std::make_pair(events, weights);
}


std::pair<std::vector<float>, std::vector<int> >
sample_pdf(TH1* hist, long int nsamples, long int maxsamples) {
  std::vector<float> events;
  std::vector<int> weights;

  double histint = hist->Integral();
  int totalbins = hist->GetNbinsX() * hist->GetNbinsY() * hist->GetNbinsZ();
  double maxperbin = maxsamples / (1.0 * totalbins);

  std::vector<int> histweights(totalbins);
  bool weighted = false;
  if (nsamples > maxsamples) {
    weighted = true;
    std::cout << "Creating weighted dataset wherever bincount > "
              << maxperbin << std::endl;

    for (int i=0; i<hist->GetNbinsX(); i++) {
      for (int j=0; j<hist->GetNbinsY(); j++) {
        for (int k=0; k<hist->GetNbinsZ(); k++) {
          double height = hist->GetBinContent(i+1, j+1, k+1);
          int idx = \
            k + j*hist->GetNbinsZ() + i*hist->GetNbinsY()*hist->GetNbinsZ();
          if (nsamples * (height / histint) > maxperbin) {
            int weight = ((nsamples * (height / histint)) / maxperbin) + 1;
            hist->SetBinContent(i+1, j+1, k+1, height/weight);
            histweights[idx] = weight;
          }
          else {
            histweights[idx] = 1.0;
          }
        }
      }
    }
  }

  int allocatesize = nsamples;
  if (weighted) {
    double totalweight = 0;
    for (int i=0; i<hist->GetNbinsX(); i++) {
      for (int j=0; j<hist->GetNbinsY(); j++) {
        for (int k=0; k<hist->GetNbinsZ(); k++) {
          int idx = \
            k + j*hist->GetNbinsZ() + i*hist->GetNbinsY()*hist->GetNbinsZ();
          totalweight += \
            hist->GetBinContent(i+1, j+1, k+1) * histweights[idx];
        }
      }
    }
    allocatesize = \
      static_cast<int>(nsamples * (1.5 * hist->Integral() / histint));
  }

  if (hist->IsA() == TH1D::Class()) {
    events.reserve(allocatesize);
    weights.reserve(allocatesize);
    TH1D* ht = dynamic_cast<TH1D*>(hist);

    double obs;
    long int j;
    for (j=0; j<nsamples; j++) {
      obs = ht->GetRandom();

      events.push_back(obs);

      if (weighted) {
        int ibin = (hist->GetXaxis()->FindBin(obs)-1);
        weights.push_back(histweights[ibin]);
        j += histweights[ibin]-1;
      }
      else {
        weights.push_back(1);
      }
    }

    if (j > nsamples) {
      weights.back() -= j-nsamples;
    }
  }
  else if (hist->IsA() == TH2D::Class()) {
    events.reserve(allocatesize * 2);
    weights.reserve(allocatesize);
    TH2D* ht = dynamic_cast<TH2D*>(hist);

    double obs0;
    double obs1;
    long int j;
    for (j=0; j<nsamples; j++) {
      ht->GetRandom2(obs0, obs1);

      events.push_back(obs0);
      events.push_back(obs1);

      if (weighted) {
        int ibin = ((hist->GetYaxis()->FindBin(obs1)-1) +
                    (hist->GetXaxis()->FindBin(obs0)-1) * hist->GetNbinsY());
        weights.push_back(histweights[ibin]);
        j += histweights[ibin]-1;
      }
      else {
        weights.push_back(1);
      }
    }
    if (j > nsamples) {
      weights.back() -= j-nsamples;
    }
  }
  else if (hist->IsA() == TH3D::Class()) {
    events.reserve(allocatesize*3);
    weights.reserve(allocatesize);
    TH3D* ht = dynamic_cast<TH3D*>(hist);

    double obs0;
    double obs1;
    double obs2;
    long int j;
    for (j=0; j<nsamples; j++) {
      ht->GetRandom3(obs0, obs1, obs2);

      events.push_back(obs0);
      events.push_back(obs1);
      events.push_back(obs2);

      if (weighted) {
        int ibin = ((hist->GetZaxis()->FindBin(obs2)-1) + 
                    (hist->GetYaxis()->FindBin(obs1)-1) * hist->GetNbinsZ() +
                    (hist->GetXaxis()->FindBin(obs0)-1) *
                     hist->GetNbinsZ()*hist->GetNbinsY());
        weights.push_back(histweights[ibin]);
        j += histweights[ibin]-1;
      }
      else {
        weights.push_back(1);
      }
    }
    if (j > nsamples) {
      weights.back() -= j-nsamples;
    }
  }
  else {
    std::cerr << "make_fake_dataset: Unknown histogram class: "
      << hist->ClassName() << std::endl;
    assert(false);
  }

  std::cout << "make_fake_dataset: Made " << weights.size()
            << ", weighted to " << nsamples << std::endl;
  return std::make_pair(events,weights);
}


unsigned nint(float nexpected) {
  return TMath::Nint(nexpected);
}

