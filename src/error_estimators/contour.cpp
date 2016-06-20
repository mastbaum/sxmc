#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <TEventList.h>
#include <TMath.h>
#include <TNtuple.h>
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
  delete this->contour_points;
}


Interval Contour::get_interval(std::string name) {
  Interval interval;
  interval.cl = this->cl;
  interval.one_sided = false;
  interval.point_estimate = -999;
  interval.coverage = -999;

  TNtuple* ls = this->contour_points;
  assert(ls);

  // Get the parameter at maximum-likelihood point
  float lmin = ls->GetMinimum("likelihood");
  float dnll = 0.13;
  long npts = 0;
  TEventList* el = NULL;
  do {
    std::ostringstream sel;
    sel << "likelihood+" << -lmin << "<" << dnll;
    ls->Draw(">>_el", sel.str().c_str());
    el = (TEventList*) gDirectory->Get("_el");
    assert(el);
    npts = el->GetN();
    dnll *= 5;
  } while (npts < 1);

  ls->SetEventList(el);
  interval.point_estimate = (ls->GetMinimum(name.c_str()) +
                             ls->GetMaximum(name.c_str())) / 2;

  ls->SetEventList(NULL);
  delete el;

  // Get the parameter extents within the contour
  interval.lower = ls->GetMinimum(name.c_str());
  interval.upper = ls->GetMaximum(name.c_str());

  ls->SetEventList(NULL);

  return interval;
}

  }  // namespace errors
}  // namespace sxmc

