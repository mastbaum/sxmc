#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <assert.h>
#include <TNtuple.h>
#include <TEnv.h>
#include <TH1F.h>
#include <TDirectory.h>

#include <sxmc/likelihood.h>
#include <sxmc/errors.h>
#include <sxmc/utils.h>

LikelihoodSpace::LikelihoodSpace(TNtuple* _samples) {
  this->samples = _samples;
  this->ml_params = extract_best_fit();
}


LikelihoodSpace::~LikelihoodSpace() {
  samples->Delete();
}


std::map<std::string, Interval> LikelihoodSpace::get_best_fit() {
  return ml_params;
}


void LikelihoodSpace::print_best_fit() {
  std::cout << "-- Best fit --" << std::endl;
  std::map<std::string, Interval>::iterator it;
  for (it=this->ml_params.begin(); it!=ml_params.end(); ++it) {
    if (it->first == "likelihood") {
      continue;
    }

    Interval interval = it->second;
    float lower_error = interval.point_estimate - interval.lower;
    float upper_error = interval.upper - interval.point_estimate;

    std::cout << " " << it->first << ": " << interval.point_estimate;
    if (interval.one_sided) {
      std::cout << " <" << interval.upper << " (" << 100 * interval.cl
                << "\% CL)" << std::endl;
    }
    else {
      std::cout << " -" << lower_error << " +" << upper_error << std::endl;
    }
  }
}


void LikelihoodSpace::print_correlations() {
  std::cout << "-- Correlation matrix --" << std::endl;
  std::vector<float> correlations = get_correlation_matrix(this->samples);

  std::vector<std::string> names;
  for (int i=0; i<this->samples->GetListOfBranches()->GetEntries(); i++) {
    std::string name = this->samples->GetListOfBranches()->At(i)->GetName();
    if (name == "likelihood") {
      continue;
    }
    names.push_back(name);
  }

  for(size_t i=0; i<names.size(); i++) {
    std::cout << std::setw(20) << names[i] << " ";
    for (size_t j=0; j<names.size(); j++) {
      std::cout << std::setiosflags(std::ios::fixed)
                << std::setprecision(3) << std::setw(8)
                << correlations.at(j + i * names.size());
    }
    std::cout << std::resetiosflags(std::ios::fixed) << std::endl;
  }
}


TH1F* LikelihoodSpace::get_projection(std::string name) {
  int default_nbins = 100;
  gEnv->GetValue("Hist.Binning.1D.x", default_nbins);
  gEnv->SetValue("Hist.Binning.1D.x", 10000);
  this->samples->Draw((name + ">>_hp").c_str(), "", "goff");
  gEnv->SetValue("Hist.Binning.1D.x", default_nbins);
  TH1F* hp = dynamic_cast<TH1F*>(gDirectory->FindObject("_hp"));
  assert(hp);
  hp->SetDirectory(NULL);

  return hp;
}


std::map<std::string, Interval>
LikelihoodSpace::extract_best_fit(ErrorType error_type) {
  // Get list of branch names
  std::vector<std::string> names;
  for (int i=0; i<this->samples->GetListOfBranches()->GetEntries(); i++) {
    names.push_back(this->samples->GetListOfBranches()->At(i)->GetName());
  }

  // Extract likelihood-maximizing parameters
  float* params_branch = new float[names.size()];
  for (size_t j=0; j<names.size(); j++) {
    this->samples->SetBranchAddress(names[j].c_str(), &params_branch[j]);
  }

  float ml_branch;
  this->samples->SetBranchAddress("likelihood", &ml_branch);

  float* params = new float[names.size()];
  float ml = 1e9;
  for (int j=0; j<this->samples->GetEntries(); j++) {
    this->samples->GetEntry(j);
    if (ml_branch < ml) {
      ml = ml_branch;
      for (size_t k=0; k<names.size(); k++) {
        params[k] = params_branch[k];
      }
    }
  }

  this->samples->ResetBranchAddresses();

  // Extract errors
  ErrorEstimator* error = NULL;
  if (error_type == ERROR_PROJECTION) {
    error = new ProjectionError(this);
  }
  else if (error_type == ERROR_CONTOUR) {
    error = new ContourError(this);
  }
  else {
    std::cerr << "LikelihoodSpace::extract_best_fit: Unknown error type"
              << std::endl;
    throw(5);
  }

  std::map<std::string, Interval> best_fit;
  for (size_t i=0; i<names.size(); i++) {
    best_fit[names[i]] = error->get_interval(names[i], params[i]);
  }

  delete error;
  delete[] params;

  return best_fit;
}

