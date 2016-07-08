#include <cassert>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <TEnv.h>
#include <TDirectory.h>
#include <TH1F.h>
#include <TNtuple.h>

#include <sxmc/contour.h>
#include <sxmc/likelihood.h>
#include <sxmc/error_estimator.h>
#include <sxmc/projection.h>
#include <sxmc/utils.h>

LikelihoodSpace::LikelihoodSpace(TNtuple* _samples, float _cl,
                                 ErrorType _error_type, bool _own_samples)
    : error_type(_error_type), cl(_cl), samples(_samples),
      own_samples(_own_samples) {
  assert(this->samples);

  // Store the ML parameters and NLL value for these samples
  this->ml_params = extract_best_fit(this->nll, this->cl, this->error_type);
}


LikelihoodSpace::~LikelihoodSpace() {
  if (own_samples) {
    delete samples;
  }
}


void LikelihoodSpace::print_best_fit() const {
  std::cout << "-- Best fit --" << std::endl;
  std::map<std::string, Interval>::const_iterator it;
  for (it=this->ml_params.begin(); it!=ml_params.end(); ++it) {
    if (it->first == "likelihood") {
      continue;
    }
    std::cout << " " << it->first << ": " << it->second.str() << std::endl;
  }

  std::cout << " NLL: " << this->nll << std::endl;
}


void LikelihoodSpace::print_correlations() {
  std::cout << "-- Correlation matrix --" << std::endl;
  std::vector<float> correlations = get_correlation_matrix(this->samples);

  std::vector<std::string> names;
  int maxlen = 0;
  for (int i=0; i<this->samples->GetListOfBranches()->GetEntries(); i++) {
    std::string name = this->samples->GetListOfBranches()->At(i)->GetName();
    if (name == "likelihood") {
      continue;
    }
    names.push_back(name);
    maxlen = std::max(maxlen, (int) name.length());
  }

  for(size_t i=0; i<names.size(); i++) {
    std::cout << std::setw(maxlen) << names[i] << " ";
    for (size_t j=0; j<names.size(); j++) {
      std::cout << std::setiosflags(std::ios::fixed)
                << std::setprecision(3) << std::setw(8)
                << correlations.at(j + i * names.size());
    }
    std::cout << std::resetiosflags(std::ios::fixed) << std::endl;
  }
}


TH1F* LikelihoodSpace::get_projection(std::string name) {
  //double default_nbins = 1;
  //gEnv->GetValue("Hist.Binning.1D.x", default_nbins);
  //gEnv->SetValue("Hist.Binning.1D.x", 5000.0);
  this->samples->Draw((name + ">>_hp").c_str(), "", "");
  TH1F* hp = dynamic_cast<TH1F*>(gDirectory->FindObject("_hp"));
  assert(hp);
  hp->SetDirectory(NULL);
  //gEnv->SetValue("Hist.Binning.1D.x", default_nbins);

  return hp;
}


TNtuple* LikelihoodSpace::get_contour(float delta) {
  float lmin = this->samples->GetMinimum("likelihood");

  std::ostringstream sel;
  sel << "likelihood+" << -lmin << "<" << delta;

  TNtuple* contour_points = \
    dynamic_cast<TNtuple*>(this->samples->CopyTree(sel.str().c_str()));

  assert(contour_points && contour_points->GetEntries() > 0);

  return contour_points;
}


std::map<std::string, Interval>
LikelihoodSpace::extract_best_fit(float& lm, float cl, ErrorType error_type) {
  // Choose error calculation method
  ErrorEstimator* error = NULL;
  if (error_type == ERROR_PROJECTION) {
    error = new sxmc::errors::Projection(this, cl);
  }
  else if (error_type == ERROR_CONTOUR) {
    error = new sxmc::errors::Contour(this, cl);
  }
  else {
    std::cerr << "LikelihoodSpace::extract_best_fit: Unknown error type"
              << std::endl;
    throw(5);
  }

  // Extract 1D parameter intervals
  std::map<std::string, Interval> best_fit;

  std::vector<std::string> names;
  for (int i=0; i<this->samples->GetListOfBranches()->GetEntries(); i++) {
    std::string name = this->samples->GetListOfBranches()->At(i)->GetName();
    if (name == "likelihood") {
      continue;
    }
    best_fit[name] = error->get_interval(name);
  }

  delete error;

  lm = this->samples->GetMinimum("likelihood");

  return best_fit;
}

