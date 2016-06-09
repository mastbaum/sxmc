#include <cassert>
#include <iostream>
#include <string>
#include <TF1.h>
#include <TH1F.h>

#include <sxmc/error_estimator.h>
#include <sxmc/likelihood.h>
#include <sxmc/projection.h>

namespace sxmc {
  namespace errors {

Interval Projection::get_interval(std::string name) { 
  Interval interval;
  interval.cl = this->cl;

  // Project out the dimension of interest
  TH1F* hp = lspace->get_projection(name);
  assert(hp && hp->Integral(0, -1) > 0);

  hp->Fit("gaus", "q goff");
  double mu = hp->GetFunction("gaus")->GetParameter(1);
  double sigma = hp->GetFunction("gaus")->GetParameter(2);

  int imax = hp->FindBin(mu);
  interval.point_estimate = mu;

  if (imax < 1) {
    imax = 1;
    interval.point_estimate = hp->GetBinLowEdge(1);
  }

  int ilo = 1;
  int ihi = 0;
  double p = 0.0;

  // Is there enough probability to the west to go central?
  if (hp->Integral(0, imax) / hp->Integral(0, -1) < this->cl / 2) {
    interval.one_sided = true;
    for (int i=0; i<hp->GetNbinsX()+1; i++) {
      p = hp->Integral(0, i) / hp->Integral(0, -1);
      if (p >= this->cl) {
        ihi = i;
        break;
      }
    }
  }
  else {
    interval.one_sided = false;

    // Find the lower limit
    for (int i=imax; i>0; i--) {
      p = hp->Integral(i, imax) / hp->Integral(0, -1);
      if (p >= this->cl / 2) {
        ilo = i;
        break;
      }
    }

    // Find the upper limit
    for (int i=imax+1;i<hp->GetNbinsX()+1; i++) {
      p = hp->Integral(imax+1, i) / hp->Integral(0, -1);
      if (p >= this->cl / 2) {
        ihi = i;
        break;
      }
    }
  }

  interval.coverage = hp->Integral(ilo, ihi) / hp->Integral(0, -1);
  interval.lower = hp->GetBinLowEdge(ilo);
  interval.upper = hp->GetBinLowEdge(ihi) + hp->GetBinWidth(ihi);

  hp->Delete();

  return interval;
}

  }  // namespace errors
}  // namespace sxmc

