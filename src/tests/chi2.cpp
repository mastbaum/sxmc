#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMath.h>
#include <TNtuple.h>

#include <sxmc/likelihood.h>
#include <sxmc/test.h>
#include <sxmc/chi2.h>

#include <sxmc/signal.h>
#include <sxmc/observable.h>
#include <sxmc/systematic.h>
#include <sxmc/source.h>

namespace sxmc {
  namespace tests {

TestResult Chi2::operator()(const LikelihoodSpace* ls) {
  std::vector<TH1*> totals(datasets->size(), NULL);

  std::map<std::string, Interval> best_fit = ls->get_best_fit();
  std::vector<float> params(best_fit.size());

  for (size_t i=0; i<sources->size(); i++) {
    std::string name = sources->at(i).name;
    params[i] = best_fit[name].point_estimate;
  }

  size_t idx = sources->size();
  for (size_t i=0; i<systematics->size(); i++) {
    for (size_t j=0; j<systematics->at(i).npars; j++) {
      std::ostringstream oss;
      oss << systematics->at(i).name << "_" << j;
      params[idx++] = best_fit[oss.str()].point_estimate;
    }
  }

  hemi::Array<unsigned> norms_buffer(signals->size(), true);
  norms_buffer.writeOnlyHostPtr();

  hemi::Array<double> param_buffer(params.size(), true);
  for (size_t i=0; i<params.size(); i++) {
    param_buffer.writeOnlyHostPtr()[i] = params[i];
  }

  // Loop over signals to build the total fit from scaled PDFs
  for (size_t i=0; i<signals->size(); i++) {
    pdfz::EvalHist* phist = \
      dynamic_cast<pdfz::EvalHist*>(signals->at(i).histogram);

    phist->SetParameterBuffer(&param_buffer, sources->size());
    phist->SetNormalizationBuffer(&norms_buffer, i);

    phist->EvalAsync(false);
    phist->EvalFinished();

    double eff = 1.0 * norms_buffer.readOnlyHostPtr()[i] / signals->at(i).n_mc;
    double nexp = \
      signals->at(i).nexpected * eff * params[signals->at(i).source.index];

    TH1* hpdf_nd = phist->CreateHistogram();

    hpdf_nd->Sumw2();
    hpdf_nd->Scale(nexp / hpdf_nd->Integral());

    unsigned ds = signals->at(i).dataset;

    if (totals[ds] == NULL) {
      totals[ds] = hpdf_nd;
    }
    else {
      totals[ds]->Add(hpdf_nd);
    }
  }

  TestResult result = TEST_OK;

  // Loop over data and compute Chi2s
  for (std::set<unsigned>::iterator it=datasets->begin();
       it!=datasets->end(); ++it) {
    unsigned ds = *it;

    TH1* total = totals[ds];
    assert(total);

    TH1* hdata = dynamic_cast<TH1*>(total->Clone("_hdata"));
    assert(hdata);
    hdata->Reset();

    if (total->IsA() == TH1D::Class()) {
      for (size_t i=0; i<samples->size()/2; i++) {
        unsigned dds = samples->at(i * 2 + 1);
        if (dds == ds) {
          ((TH1D*) hdata)->Fill(samples->at(i * 2 + 0));
        }
      }
    }
    else if (total->IsA() == TH2D::Class()) {
      for (size_t i=0; i<samples->size()/3; i++) {
        unsigned dds = samples->at(i * 3 + 2);
        if (dds == ds) {
          ((TH2D*) hdata)->Fill(samples->at(i * 3 + 0), samples->at(i * 3 + 1));
        }
      }
    }
    else if (total->IsA() == TH3D::Class()) {
      for (size_t i=0; i<samples->size()/4; i++) {
        unsigned dds = samples->at(i * 4 + 3);
        if (dds == ds) {
          ((TH3D*) hdata)->Fill(samples->at(i * 4 + 0), samples->at(i * 4 + 1), samples->at(i * 4 + 2));
        }
      }
    }
    else {
      std::cerr << "tests::Chi2: Unknown histogram type" << std::endl;
      assert(false);
    }
    
    double chi2ndf = hdata->Chi2Test(total, "UW CHI2/NDF");

    std::cout << "tests::Chi2: Dataset " << ds << ": "
              << "Chi2/ndf = " << chi2ndf << std::endl;

    if (chi2ndf > p_threshold) {
      result = TEST_FAIL;
    }
  }

  std::cout << "tests::Chi2: " << status(result) << std::endl;

  return result;
}

  }  // namespace tests
}  // namespace sxmc

