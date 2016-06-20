#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TH1F.h>
#include <TMath.h>
#include <TNtuple.h>

#include <sxmc/likelihood.h>
#include <sxmc/test.h>
#include <sxmc/drift.h>

namespace sxmc {
  namespace tests {

TestResult Drift::operator()(const LikelihoodSpace* ls) {
  TNtuple* samples = const_cast<TNtuple*>(ls->get_samples());
  assert(samples && samples->GetEntries() > 2 * this->slices);

  TestResult result = TEST_OK;

  unsigned chunk = 1.0 * samples->GetEntries() / this->slices;

  // Loop over parameters
  for (int i=0; i<samples->GetListOfBranches()->GetEntries(); i++) {
    std::string name = samples->GetListOfBranches()->At(i)->GetName();

    std::vector<size_t> ns;
    std::vector<double*> points;

    // Loop over time slices
    for (unsigned i=0; i<slices; i++) {
      //std::ostringstream sname;
      //sname << "_" << name << "_" << i;

      //TH1F* h = NULL;

      // Clone the first, to preserve binning
      //if (i > 0) {
      //  h = dynamic_cast<TH1F*>(projections[0]->Clone(sname.str().c_str()));
      //  assert(h);
      //  h->Reset();
      //}

      unsigned start = 1.0 * i * chunk;
      unsigned end = 1.0 * (i + 1) * chunk;

      std::ostringstream scut;
      scut << "Entry$>=" << start << "&&Entry$<" << end;

      long n = samples->Draw(name.c_str(), //(name + ">>" + sname.str()).c_str(),
                    scut.str().c_str(), "goff");
      double* vals = samples->GetVal(0);
      double* v2 = new double[n];

      memcpy(v2, vals, n * sizeof(double));

      std::sort(v2, v2+n);

      //h = dynamic_cast<TH1F*>(gDirectory->Get(sname.str().c_str()));
      //assert(h);

      //if (i == 0) {
      //  h->Rebin(4);
      //}

      //TCanvas c1;
      //h->Draw();
      //c1.SaveAs((sname.str() + ".pdf").c_str());

      ns.push_back(n);
      points.push_back(v2);
    }

    // Compare each slice to each other, computing a Chi^2 p-value
    for (size_t i=0; i<points.size()-1; i++) {
      for (size_t j=i+1; j<points.size(); j++) {
        //TH1F* a = projections[i];
        //TH1F* b = projections[j];

        double p = TMath::KolmogorovTest(ns[i], points[i], ns[j], points[j], "");

        if (p < this->p_threshold) {
          std::cout << "tests::Drift: " << name << " "
                    << "slices " << i << "/" << j << ": "
                    << "p = " << p << ", threshold = " << this->p_threshold
                    << std::endl;
          result = TEST_FAIL;
        }

        //if (a->GetEntries() < 1000 || b->GetEntries() < 1000) {
        //  std::cout << "tests::Drift: Insufficient statistics for "
        //            << name << " slices " << i << "/" << j << std::endl;
        //  result = TEST_UNKNOWN;
        //  continue;
        //}

        //double p = a->Chi2Test(b,"UU P CHI2/NDF");
        //if (p < this->p_threshold) {
        //  std::cout << "tests::Drift: " << name << " "
        //            << "slices " << i << "/" << j << ": "
        //            << "p = " << p << ", threshold = " << this->p_threshold
        //            << std::endl;
        //  result = TEST_FAIL;
        //}
      }
    }

    for (size_t i=0; i<points.size(); i++) {
      delete[] points[i];
    }
  }

  std::cout << "tests::Drift: " << status(result) << std::endl;

  return result;
}

  }  // namespace tests
}  // namespace sxmc

