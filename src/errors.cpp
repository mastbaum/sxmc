#include <string>
#include <sstream>
#include <TH1F.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TEnv.h>
#include <TDirectory.h>
#include <assert.h>

#include <sxmc/likelihood.h>
#include <sxmc/errors.h>

std::string Interval::str() {
  float lower_error = this->point_estimate - this->lower;
  float upper_error = this->upper - this->point_estimate;

  std::stringstream ss;
  ss << this->point_estimate;
  if (this->one_sided) {
    ss << " <" << this->upper << " (" << 100 * this->cl << "\% CL)";
  }
  else {
    ss << " -" << lower_error << " +" << upper_error;
  }

  return ss.str();
}


Interval ProjectionError::get_interval(std::string name,
                                       float point_estimate) { 
  // Project out the dimension of interest
  TH1F* hp = lspace->get_projection(name);

  float integral = hp->Integral(0, hp->GetNbinsX());
  assert(integral > 0);

  Interval interval;
  interval.point_estimate = point_estimate;
  interval.cl = this->cl;

  float integrated_prob = 0;  // for coverage calculation
  float alpha = (1.0 - cl) / 2;

  if (hp->GetMean() < 2 * hp->GetRMS()) {
    interval.lower = 0;
    interval.one_sided = true;
    alpha *= 2;
  }
  else {
    // Find the lower bound
    int lower_bin = 0;
    for (int i=0; i<hp->GetNbinsX(); i++) {
      if (hp->Integral(0, i) / integral >= alpha) {
        lower_bin = i;
        break;
      }
    }
    interval.lower = hp->GetBinLowEdge(lower_bin);
    integrated_prob += hp->Integral(0, lower_bin);
    interval.one_sided = false;
  }

  // Find the upper bound
  int upper_bin = hp->GetNbinsX();
  for (int i=hp->GetNbinsX(); i>0; i--) {
    if (hp->Integral(i, hp->GetNbinsX()) / integral > alpha) {
      upper_bin = i;
      break;
    }
  }
  interval.upper = hp->GetBinLowEdge(upper_bin);
  integrated_prob += hp->Integral(upper_bin, hp->GetNbinsX());

  // Coverage
  interval.coverage = integrated_prob / integral;

  hp->Delete();

  return interval;
}


ContourError::ContourError(LikelihoodSpace* _lspace, float _cl)
    : ErrorEstimator(_lspace, _cl) {
  float delta = 0.5 * TMath::ChisquareQuantile(_cl, 1);
  contour_points = _lspace->get_contour(delta);
}


ContourError::~ContourError() {
  contour_points->Delete();
}


Interval ContourError::get_interval(std::string name, float point_estimate) {
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
  assert(hp);
  hp->SetDirectory(NULL);

  // Find extrema
  int ll_bin = hp->FindFirstBinAbove(0);
  interval.lower = hp->GetBinLowEdge(ll_bin);

  int ul_bin = hp->FindLastBinAbove(0);
  interval.upper = hp->GetBinLowEdge(ul_bin) +
                   hp->GetBinWidth(ul_bin);

  return interval;
}

/*
    // parameters for counts -> lifetime -> mass conversion
    float n_te130 = 7.46e26 * TMath::Power(radius_cut / 3500, 3);
    float Mbb = 4.03;  // IBM-2
    float Gphase = 3.69e-14;  // y^-1, using g_A = 1.269
    float m_beta = 511e3;  // 511 keV, in eV

    // calculate sensitivity by integrating over backgrounds
    float lifetime = n_te130 * live_time * 0.69315 * signal_eff / limit;
    float mass = m_beta / TMath::Sqrt(lifetime * Gphase * Mbb * Mbb);

    std::cout << "-- Sensitivity (marginalized backgrounds): "
              << confidence * 100 << "\% CL --" << std::endl;
    std::cout << " Counts: " << limit << std::endl;
    std::cout << " T_1/2:  " << lifetime << " y" << std::endl;
    std::cout << " Mass:   " << 1000 * mass << " meV" << std::endl;

    limits["counts"].Fill(limit);
    limits["lifetime"].Fill(lifetime);
    limits["mass"].Fill(1000 * mass);

    // sensitivity using minos-like contour
    float lifetime_contour = \
      n_te130 * live_time * 0.69315 * signal_eff / ul_contour;
    float mass_contour = \
      m_beta / TMath::Sqrt(lifetime_contour * Gphase * Mbb * Mbb);

    std::cout << "-- Contour Sensitivity: " << confidence * 100 << "\% CL --"
              << std::endl;
    std::cout << " Delta Chi2: " << 2 * delta << std::endl;
    std::cout << " UL Counts:  " << ul_contour << std::endl;
    std::cout << " LL Counts:  " << ll_contour << std::endl;
    std::cout << " T_1/2:      " << lifetime_contour << " y" << std::endl;
    std::cout << " Mass:       " << 1000 * mass_contour << " meV" << std::endl;

    limits["contour_counts"].Fill(ul_contour);
    limits["contour_lifetime"].Fill(lifetime_contour);
    limits["contour_mass"].Fill(1000 * mass_contour);
    if (0 >= ll_contour && 0 <= ul_contour) {
      std::cout << "covered" << std::endl;
      hcontcoverage.Fill(nevents);
    }

    // sensitivity using feldman-cousins
    if (nevents < 25) {
      TFeldmanCousins tfc(confidence);
      tfc.SetMuMax(2500);
      float fc_limit = tfc.CalculateUpperLimit(nevents, nexpected);
      float fc_llimit = tfc.CalculateLowerLimit(nevents, nexpected);
      float fc_lifetime = n_te130 * live_time * 0.69315 * signal_eff / fc_limit;
      float fc_mass = m_beta / TMath::Sqrt(fc_lifetime * Gphase * Mbb * Mbb);
      std::cout << "-- Feldman-Cousins Sensitivity: "
                << 100 * confidence << "\% CL --" << std::endl;
      std::cout << " UL Counts: " << fc_limit << std::endl;
      std::cout << " LL Counts: " << fc_llimit << std::endl;
      std::cout << " T_1/2:     " << fc_lifetime << " y" << std::endl;
      std::cout << " Mass:      " << 1000 * fc_mass << " meV" << std::endl;
  
      limits["fc_counts"].Fill(fc_limit);
      limits["fc_lifetime"].Fill(fc_lifetime);
      limits["fc_mass"].Fill(1000 * fc_mass);
      if (0 >= fc_llimit && 0 <= fc_limit) {
        std::cout << "covered" << std::endl;
        hfccoverage.Fill(nevents);
      }
    }
    else {
      std::cout << "Skipping Feldman-Cousins calculation" << std::endl;
    }
*/

