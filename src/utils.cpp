#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1.h>
#include <TDirectory.h>
#include <TNtuple.h>
#include <TMath.h>

#include <sxmc/utils.h>

unsigned nint(float nexpected) {
  return TMath::Nint(nexpected);
}


float get_ntuple_entry(TNtuple* nt, int i, std::string field) {
  float v;
  nt->SetBranchAddress(field.c_str(), &v);
  assert(i < nt->GetEntries());
  nt->GetEvent(i);
  nt->ResetBranchAddresses();
  return v;
}


std::vector<float> get_correlation_matrix(TNtuple* nt) {
  int nentries = nt->GetEntries();

  // Get list of branch names
  std::vector<std::string> names;
  for (int i=0; i<nt->GetListOfBranches()->GetEntries(); i++) {
    std::string name = nt->GetListOfBranches()->At(i)->GetName();
    if (name == "likelihood") {
      continue;
    }
    names.push_back(name);
  }

  std::vector<float> matrix(names.size() * names.size());

  // Convert the ntuple to a vector, calculating means as we go
  std::vector<float> table(names.size() * nentries);
  std::vector<float> means(names.size(), 0);
  for (int i=0; i<nentries; i++) {
    for (size_t j=0; j<names.size(); j++) {
      float v = get_ntuple_entry(nt, i, names.at(j));
      table.at(j + i * names.size()) = v;
      means.at(j) += v;
    }
  }

  // Sums to means
  for (size_t i=0; i<names.size(); i++) {
    means.at(i) /= nentries;
  }

  // Compute correlations
  for (size_t i=0; i<names.size(); i++) {
    for (size_t j=i; j<names.size(); j++) {
      float t = 0;
      float dx2 = 0;
      float dy2 = 0;
      for (int k=0; k<nentries; k++) {
        float x1 = table.at(i + k * names.size()) - means.at(i);
        float x2 = table.at(j + k * names.size()) - means.at(j);
        t += x1 * x2;
        dx2 += x1 * x1;
        dy2 += x2 * x2;
      }
      matrix.at(i * names.size() + j) = t / TMath::Sqrt(dx2 * dy2);
    }
  }

  return matrix;
}

