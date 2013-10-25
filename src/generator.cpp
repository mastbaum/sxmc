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

std::vector<float> sample_pdf(TH1* hist, int nsamples)
{
  size_t obs_id = 0;
  std::vector<float> events;
  if (hist->IsA() == TH1D::Class()) {
    events.resize(nsamples);
    TH1D* ht = dynamic_cast<TH1D*>(hist);

    double obs;
    for (int j=0; j<nsamples; j++) {
      obs = ht->GetRandom();

      events[obs_id++] = obs;
    }
  }
  else if (hist->IsA() == TH2D::Class()) {
    events.resize(nsamples*2);
    TH2D* ht = dynamic_cast<TH2D*>(hist);

    double obs0;
    double obs1;
    for (int j=0; j<nsamples; j++) {
      ht->GetRandom2(obs0, obs1);

      events[obs_id++] = obs0;
      events[obs_id++] = obs1;
    }
  }
  else if (hist->IsA() == TH3D::Class()) {
    events.resize(nsamples*3);
    TH3D* ht = dynamic_cast<TH3D*>(hist);

    double obs0;
    double obs1;
    double obs2;
    for (int j=0; j<nsamples; j++) {
      ht->GetRandom3(obs0, obs1, obs2);

      events[obs_id++] = obs0;
      events[obs_id++] = obs1;
      events[obs_id++] = obs2;
    }
  }
  else {
    std::cerr << "make_fake_dataset: Unknown histogram class: "
      << hist->ClassName() << std::endl;
    assert(false);
  }
  return events;
}

std::vector<float> make_fake_dataset(std::vector<Signal>& signals,
                                     std::vector<Systematic>& systematics,
                                     std::vector<Observable>& observables,
                                     std::vector<float> params, bool poisson) {
  std::cout << "make_fake_dataset: Generating dataset..." << std::endl;

  std::vector<float> syst_vals;
  for (size_t i=signals.size();i<signals.size()+systematics.size();i++){
    syst_vals.push_back(params[i]); 
  }
  std::vector<float> upper;
  std::vector<float> lower;
  for (size_t i=0;i<observables.size();i++){
    upper.push_back(observables[i].upper);
    lower.push_back(observables[i].lower);
  }

  std::vector<float> events;
  std::vector<unsigned> observed(signals.size());
  for (size_t i=0;i<signals.size();i++){
    observed[i] = \
      dynamic_cast<pdfz::EvalHist*>(signals[i].histogram)->RandomSample(events,params[i],
          syst_vals,upper,lower,poisson);

    std::cout << "make_fake_dataset: " << signals[i].name << ": "
              << observed[i] << " events (" << signals[i].nexpected
              << " expected)" << std::endl;
  }

  return events;
}


unsigned nint(float nexpected)
{
  return TMath::Nint(nexpected);
}

