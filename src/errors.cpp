#include <string>
#include <TH1F.h>
#include <assert.h>

#include <sxmc/likelihood.h>
#include <sxmc/errors.h>

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

/*
    // find contour
    std::cout << "Generating " << 100 * confidence
              << "\% CL contour..." << std::endl;
    float delta = 0.5 * TMath::ChisquareQuantile(confidence, 1);
    TNtuple* lscontour = (TNtuple*) lspace->Clone("lscontour");
    lscontour->Reset();

    float* v = new float[params.size() + 1];
    for (int j=0; j<lspace->GetEntries(); j++) {
      lspace->GetEntry(j);
      // build a list of points inside the contour
      if (ml_branch < ml + delta) {
        for (size_t k=0; k<params.size(); k++) {
          v[k] = params_branch[k];
        }
        v[params.size()] = ml_branch;
        lscontour->Fill(v);
      }
    }
    delete[] v;
    
    TH1F hcproj("hcproj", "hcproj", 20000, 0, 500);
    lscontour->Draw((signal_name + ">>hcproj").c_str(), "", "goff");
    int ul_bin = hcproj.FindLastBinAbove(0);
    int ll_bin = hcproj.FindFirstBinAbove(0);
    double ul_contour = hcproj.GetBinLowEdge(ul_bin) +
                        hcproj.GetBinWidth(ul_bin);
    double ll_contour = hcproj.GetBinLowEdge(ll_bin);

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

