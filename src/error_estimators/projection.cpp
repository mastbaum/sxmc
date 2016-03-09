#include <assert.h>
#include <iostream>
#include <string>
#include <TH1F.h>

#include <sxmc/error_estimator.h>
#include <sxmc/likelihood.h>
#include <sxmc/projection.h>

namespace sxmc {
  namespace errors {

Interval Projection::get_interval(std::string name,
                                  float point_estimate) { 
  // Project out the dimension of interest
  TH1F* hp = lspace->get_projection(name);

  float integral = hp->Integral(0, hp->GetNbinsX() + 1);
  assert(integral > 0);

  Interval interval;
  interval.point_estimate = point_estimate;
  interval.cl = this->cl;

  float integrated_prob = 0;  // for coverage calculation
  float alpha = (1.0 - cl) / 2;

  int lower_bin = 0;
  if (hp->GetMean() < 2 * hp->GetRMS()) {
    interval.lower = 0;
    interval.one_sided = true;
    alpha *= 2;
  }
  else {
    // Find the lower bound
    for (int i=0; i<hp->GetNbinsX(); i++) {
      if (hp->Integral(0, i) / integral >= alpha) {
        lower_bin = i;
        break;
      }
    }
    interval.lower = hp->GetBinLowEdge(lower_bin);
    interval.one_sided = false;
  }

  // Find the upper bound
  int upper_bin = hp->GetNbinsX() + 1;
  for (int i=hp->GetNbinsX()+1; i>0; i--) {
    if (hp->Integral(i, hp->GetNbinsX() + 1) / integral > alpha) {
      upper_bin = i;
      break;
    }
  }
  interval.upper = hp->GetBinLowEdge(upper_bin);

  // Coverage
  interval.coverage = hp->Integral(lower_bin, upper_bin) / integral;

  hp->Delete();

  return interval;
}

  }  // namespace errors
}  // namespace sxmc

