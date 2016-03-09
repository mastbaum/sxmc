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

#include <sxmc/contour.h>
#include <sxmc/likelihood.h>
#include <sxmc/error_estimator.h>
#include <sxmc/projection.h>
#include <sxmc/utils.h>

LikelihoodSpace::LikelihoodSpace(TNtuple* _samples, float cl) {
  this->samples = _samples;
  this->ml_params = extract_best_fit(this->ml, cl);
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

    std::cout << " " << it->first << ": " << it->second.str() << std::endl;
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


TNtuple* LikelihoodSpace::get_contour(float delta) {
  TNtuple* contour = (TNtuple*) this->samples->Clone("lscontour");
  contour->Reset();

  // Get list of branch names
  std::vector<std::string> names;
  for (int i=0; i<this->samples->GetListOfBranches()->GetEntries(); i++) {
    std::string name = this->samples->GetListOfBranches()->At(i)->GetName();
    if (name == "likelihood") {
      continue;
    }
    names.push_back(name);
  }

  float* params_branch = new float[names.size()];
  for (size_t i=0; i<names.size(); i++) {
    this->samples->SetBranchAddress(names[i].c_str(), &params_branch[i]);
  }

  float ml_branch;
  this->samples->SetBranchAddress("likelihood", &ml_branch);

  // Build a new TNtuple with samples inside the contour
  float* v = new float[names.size() + 1];
  for (int i=0; i<this->samples->GetEntries(); i++) {
    this->samples->GetEntry(i);
    if (ml_branch < this->ml + delta) {
      for (size_t j=0; j<names.size(); j++) {
        v[j] = params_branch[j];
      }
      v[names.size()] = ml_branch;
      contour->Fill(v);
    }
  }

  this->samples->ResetBranchAddresses();
  delete[] v;

  return contour;
}


std::map<std::string, Interval>
LikelihoodSpace::extract_best_fit(float& ml, float cl, ErrorType error_type) {
  // Get list of branch names
  std::vector<std::string> names;
  for (int i=0; i<this->samples->GetListOfBranches()->GetEntries(); i++) {
    std::string name = this->samples->GetListOfBranches()->At(i)->GetName();
    if (name == "likelihood") {
      continue;
    }
    names.push_back(name);
  }

  // Extract likelihood-maximizing parameters
  float* params_branch = new float[names.size()];
  for (size_t j=0; j<names.size(); j++) {
    this->samples->SetBranchAddress(names[j].c_str(), &params_branch[j]);
  }

  float ml_branch;
  this->samples->SetBranchAddress("likelihood", &ml_branch);

  float* params = new float[names.size()];
  ml = 1e9;
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

  std::map<std::string, Interval> best_fit;
  for (size_t i=0; i<names.size(); i++) {
    best_fit[names[i]] = error->get_interval(names[i], params[i]);
  }

  delete error;
  delete[] params;

  return best_fit;
}

