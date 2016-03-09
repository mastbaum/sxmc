#include <iostream>
#include <string>
#include <TH1F.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TEnv.h>
#include <TDirectory.h>
#include <assert.h>

#include <sxmc/contour.h>
#include <sxmc/error_estimator.h>
#include <sxmc/likelihood.h>

namespace sxmc {
  namespace errors {

Contour::Contour(LikelihoodSpace* _lspace, float _cl)
    : ErrorEstimator(_lspace, _cl) {
  float delta = 0.5 * TMath::ChisquareQuantile(_cl, 1);
  contour_points = _lspace->get_contour(delta);
}


Contour::~Contour() {
  contour_points->Delete();
}


Interval Contour::get_interval(std::string name, float point_estimate) {
  Interval interval;
  interval.cl = this->cl;
  interval.one_sided = false;
  interval.point_estimate = point_estimate;
  interval.coverage = -1;

  // Plot projection of the requested parameter from the reduced space
  int default_nbins = 100;
  gEnv->GetValue("Hist.Binning.1D.x", default_nbins);
  gEnv->SetValue("Hist.Binning.1D.x", 20000);
  this->contour_points->Draw((name + ">>_hp").c_str(), "", "goff");
  gEnv->SetValue("Hist.Binning.1D.x", default_nbins);
  TH1F* hp = dynamic_cast<TH1F*>(gDirectory->FindObject("_hp"));

  if (hp) {
    hp->SetDirectory(NULL);

    // Find extrema
    int ll_bin = hp->FindFirstBinAbove(0);
    interval.lower = hp->GetBinLowEdge(ll_bin);

    int ul_bin = hp->FindLastBinAbove(0);
    interval.upper = hp->GetBinLowEdge(ul_bin) +
                     hp->GetBinWidth(ul_bin);
  }
  else {
    std::cerr << "Contour::get_interval: Error getting projection"
              << std::endl;
  }

  delete hp;
  return interval;
}

  }  // namespace errors
}  // namespace sxmc

