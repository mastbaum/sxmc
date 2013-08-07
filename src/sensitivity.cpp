/**
 * \file sensitivity.cpp
 * \brief Calculate sensitivity with an MCMC
 *
 * An ensemble of fake experiments is generated and each fit with an MCMC, with
 * parameters defined in the JSON configuration file. Limits are calculated for
 * each.
 */

#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TRandom.h>
#include <TRandom2.h>
#include <TMath.h>
#include <TFeldmanCousins.h>
#include "config.h"
#include "generator.h"
#include "mcmc.h"
#include "utils.h"

/** A ROOT color palette in which sequential colors look okay. */
static const int color[28] = {kRed,      kGreen,    kBlue,      kMagenta,
                              kCyan,     kYellow,   kOrange,    kViolet+2,
                              kRed+2,    kGreen+2,  kBlue+2,    kMagenta+2,
                              kCyan+2,   kYellow+2, kOrange+2,  kRed-7,
                              kGreen-7,  kBlue-7,   kMagenta-7, kCyan-7,
                              kYellow-7, kOrange-7, kRed-6,     kAzure+1,
                              kTeal+1,   kSpring-9, kAzure-9};


/**
 * Find an upper limit
 *
 * Locate the x value that a fraction 1-CL/2 of the distribution falls above.
 * This may over-cover due to finite histogram bin widths.
 *
 * \param h 1D histogram to calculate limit on
 * \param cl Confidence level
 * \returns The nearest bin edge covering at least the desired confidence
 */
double upper_limit(TH1F* h, double cl=0.682) {
  double tail = (1.0 - cl) / 2;
  double integral = 0;
  int thisbin = h->GetNbinsX();
  while (integral < tail * h->Integral()) {
    integral = h->Integral(thisbin--, h->GetNbinsX()+1);
  }
  return h->GetBinLowEdge(thisbin + 1);
}


/**
 * Find a confidence interval
 *
 * Locate the boundaries of a confidence interval, optionally bounded at zero.
 * Typically, we integrate outward from the central value such that a
 * fraction CL/2 falls on either side. In the bounded case, the interval is
 * shifted upward so that it all falls above 0, and the upper boundary is
 * reported as a limit, in the spirit of Feldman & Cousins.
 *
 * This may over-cover due to finite histogram bin widths.
 *
 * \param h 1D histogram to calculate limit on
 * \param mu Central value
 * \param lower Lower boundary of the interval
 * \param upper Upper boundary of the interval
 * \param limit True if the upper boundary is a limit
 * \param bounded Apply a boundary at zero
 * \param cl Confidence level
 */
void find_interval(TH1F* h, float mu, float& lower, float& upper, bool& limit,
                   bool bounded=false, double cl=0.682) {
  int central_bin = h->FindBin(mu);
  float total_integral = h->Integral();
  limit = false;

  // lower boundary
  int lower_bin = central_bin;
  float integral = 0;
  while (integral < cl/2 * total_integral) {
    integral = h->Integral(lower_bin--, central_bin);

    // crashed into zero?
    if ((h->GetBinLowEdge(lower_bin) <= 0 && bounded) || lower_bin == 0) {
      lower_bin = 1;
      limit = true;
      break;
    }
  }

  // upper boundary, starting from lower boundary
  int upper_bin = lower_bin;
  integral = 0;
  while (integral < cl * total_integral) {
    integral = h->Integral(lower_bin, upper_bin++);

    // we've run out of pdf... shift the interval down
    if (upper_bin > h->GetNbinsX() + 1) {
      upper_bin--;
      lower_bin--;
    }
  }

  lower = h->GetBinLowEdge(lower_bin);
  upper = h->GetBinLowEdge(upper_bin);
}


/**
 * Run an ensemble of independent fake experiments
 *
 * Run experiments, fit with mcmc, and tabulate background-fluctuation
 * sensitivity for each, creating a histogram of limits. The estimated
 * sensitivity is the median of the limits of the ensemble.
 *
 * \param signals List of Signals defining PDFs, rates, etc.
 * \param systematics List of Systematics applied to PDFs
 * \param observables List of Observables common to PDFs
 * \param cuts List of cuts applied to data set
 * \param steps Number of MCMC random walk steps to take
 * \param burnin_fraction Fraction of initial MCMC steps to throw out
 * \param signal_name The name of the Signal that is the signal
 * \param signal_eff Efficiency for detecting signal
 * \param confidence Desired confidence level for limits
 * \param nexperiments Number of fake experiments to run
 * \param debug_mode If true, accept and save all steps
 * \returns A (name, histogram limits from the experiments) map
 */
std::map<std::string, TH1F> ensemble(std::vector<Signal>& signals,
                                     std::vector<Systematic>& systematics,
                                     std::vector<Observable>& observables,
                                     std::vector<Observable>& cuts,
                                     unsigned steps, float burnin_fraction,
                                     std::string signal_name, float signal_eff,
                                     float confidence, unsigned nexperiments,
                                     float live_time, const bool debug_mode) {
  std::map<std::string, TH1F> limits;
  limits["counts"] = TH1F("limits_c", ";Counts in fit range;Fraction",
                          10000, 0, 500);
  limits["lifetime"] = TH1F("limits_l", ";T_{1/2} limit (y);Fraction",
                            10000, 1e24, 1e27);
  limits["mass"] = TH1F("limits_m", ";m_{#beta#beta} limit (meV);Fraction",
                        1000, 0, 1000);
  limits["fc_counts"] = TH1F("fc_limits_c", ";Counts in fit range;Fraction",
                             1000, 0, 500);
  limits["fc_lifetime"] = TH1F("fc_limits_l", ";T_{1/2} limit (y);Fraction",
                               10000, 1e24, 1e27);
  limits["fc_mass"] = TH1F("fc_limits_m",
                           ";m_{#beta#beta} limit (meV);Fraction",
                           10000, 0, 1000);
  limits["contour_counts"] = TH1F("contour_limits_c",
                                  ";Counts in fit range;Fraction",
                                  10000, 0, 500);
  limits["contour_lifetime"] = TH1F("contour_limits_l",
                                    ";T_{1/2} limit (y);Fraction",
                                    10000, 1e24, 1e27);
  limits["contour_mass"] = TH1F("contour_limits_m",
                                ";m_{#beta#beta} limit (meV);Fraction",
                                1000, 0, 1000);

  // coverage calculation
  TH1F hfccoverage("hfccoverage", ";True N;Coverage", 50, 0, 50);
  hfccoverage.Sumw2();
  TH1F hcontcoverage("hcontcoverage", ";True N;Coverage", 50, 0, 50);
  hcontcoverage.Sumw2();
  TH1F hcounts("hcounts", ";True N;Event count;", 50, 0, 50);
  hcounts.Sumw2();

  for (unsigned i=0; i<nexperiments; i++) {
    std::cout << "Experiment " << i << " / " << nexperiments << std::endl;

    // make fake data
    std::vector<float> params;
    std::vector<std::string> param_names;
    float nexpected = 0;
    for (size_t j=0; j<signals.size(); j++) {
      params.push_back(signals[j].nexpected);
      param_names.push_back(signals[j].name);
      nexpected += signals[j].nexpected;
    }
    for (size_t j=0; j<systematics.size(); j++) {
      params.push_back(systematics[j].mean);
      param_names.push_back(systematics[j].name);
    }
    std::vector<float> data = \
      make_fake_dataset(signals, systematics, observables, params, true);

    size_t nevents = (int) 1.0 * data.size() / observables.size();
    hcounts.Fill(nevents);

    // run mcmc
    MCMC mcmc(signals, systematics, observables);
    TNtuple* lspace = mcmc(data, steps, burnin_fraction, debug_mode);

    TFile f("lspace.root", "recreate");
    TNtuple* lsclone = (TNtuple*) lspace->Clone("ls");
    lsclone->Write();
    f.Close();

    // calculate signal sensitivity
    TH1F hproj("hproj", "hproj", 20000, 0, 500);
    lspace->Draw((signal_name + ">>hproj").c_str(), "", "goff");
    double limit = upper_limit(&hproj, confidence);

    float radius_cut = 0;
    for (size_t j=0; j<observables.size(); j++) {
      if (observables[j].name == "radius") {
        radius_cut = observables[j].upper;
        break;
      }
    }
    for (size_t j=0; j<cuts.size(); j++) {
      if (cuts[j].name == "radius") {
        radius_cut = cuts[j].upper;
        break;
      }
    }

    // extract likelihood-maximizing parameters
    std::cout << "Finding maximum likelihood..." << std::endl;
    float* params_branch = new float[params.size()];
    for (size_t j=0; j<params.size(); j++) {
      lspace->SetBranchAddress(param_names[j].c_str(), &params_branch[j]);
    }
    float ml_branch;
    lspace->SetBranchAddress("likelihood", &ml_branch);
    float* params_fit = new float[params.size()];
    float ml = 1e9;
    float nfit = 0;
    for (int j=0; j<lspace->GetEntries(); j++) {
      lspace->GetEntry(j);
      if (ml_branch < ml) {
        ml = ml_branch;
        nfit = 0;
        for (size_t k=0; k<params.size(); k++) {
          params_fit[k] = params_branch[k];
          if (k < signals.size()) {
            nfit += params_fit[k];
          }
        }
      }
    }
    std::cout << "Fit " << nfit << " events" << std::endl;

    // find contour
    std::cout << "Generating " << 100 * confidence
              << " contour..." << std::endl;
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

    // output best-fit parameters
    std::cout << "-- Best fit: NLL = " << ml << " --" << std::endl;
    for (size_t j=0; j<params.size(); j++) {
      float fp = params_fit[j];
      float sqrt_fp = TMath::Sqrt(fp);
      TH1F hp("hp", "hp", 10000, fp - 50.0 * sqrt_fp, fp + 50.0 * sqrt_fp);
      lspace->Draw((param_names[j] + ">>hp").c_str(), "", "goff");
      float lower_boundary;
      float upper_boundary;
      bool is_limit;
      find_interval(&hp, params_fit[j], lower_boundary, upper_boundary,
                    is_limit, (j < signals.size()), 0.682);
      float lower_error = params_fit[j] - lower_boundary;
      float upper_error = upper_boundary - params_fit[j];

      std::cout << " " << param_names[j] << ": " << params_fit[j];
      if (is_limit) {
        std::cout << " <" << upper_error << " @ 68\% CL";
      }
      else {
        std::cout << " -" << lower_error << " +" << upper_error;
      }
      std::cout << (j < signals.size() ? " (N = " : " (mean = ")
                << params[j] <<  ")";
      std::cout << std::endl;
    }

    std::cout << "-- Correlation matrix --" << std::endl;
    std::vector<float> correlations = get_correlation_matrix(lspace);
    for (size_t j=0; j<params.size(); j++) {
      std::cout << std::setw(20) << param_names[j] << " ";
      for (size_t k=0; k<params.size(); k++) {
        std::cout << std::setiosflags(std::ios::fixed)
                  << std::setprecision(3) << std::setw(8)
                  << correlations[k + j * (params.size() + 1)];
      }
      std::cout << std::resetiosflags(std::ios::fixed) << std::endl;
    }

    // plot and save this fit
    std::vector<SpectralPlot> plots_full;
    std::vector<SpectralPlot> plots_external;
    std::vector<SpectralPlot> plots_cosmogenic;
    for (size_t j=0; j<observables.size(); j++) {
      Observable* o = &observables[j];
      std::stringstream ytitle;
      ytitle << "Counts/" << (o->upper - o->lower) / o->bins << " "
             << o->units << "/" << live_time << " y";
      plots_full.push_back(SpectralPlot(2, o->lower, o->upper, 1e-2, 1e6,
                           true, "", o->title, ytitle.str().c_str()));
      plots_external.push_back(SpectralPlot(2, o->lower, o->upper, 1e-2, 1e6,
                               true, "", o->title, ytitle.str().c_str()));
      plots_cosmogenic.push_back(SpectralPlot(2, o->lower, o->upper, 1e-2, 1e6,
                                 true, "", o->title, ytitle.str().c_str()));
    }

    std::vector<TH1D*> external_total(observables.size(), NULL);
    std::vector<TH1D*> cosmogenic_total(observables.size(), NULL);
    std::vector<TH1D*> fit_total(observables.size(), NULL);

    hemi::Array<unsigned> norms_buffer(signals.size(), true);
    norms_buffer.writeOnlyHostPtr();

    hemi::Array<double> param_buffer(params.size(), true);
    for (size_t j=0; j<params.size(); j++) {
      param_buffer.writeOnlyHostPtr()[j] = params_fit[j];
    }

    for (size_t j=0; j<signals.size(); j++) {
      pdfz::EvalHist* phist = \
        dynamic_cast<pdfz::EvalHist*>(signals[j].histogram);
      phist->SetParameterBuffer(&param_buffer, signals.size());
      phist->SetNormalizationBuffer(&norms_buffer, j);

      TH1* hpdf_nd = phist->CreateHistogram();
      hpdf_nd->Scale(params_fit[j] / hpdf_nd->Integral());

      std::vector<TH1D*> hpdf(observables.size(), NULL);
      if (hpdf_nd->IsA() == TH1D::Class()) {
        hpdf[0] = dynamic_cast<TH1D*>(hpdf_nd);
      }
      else if (hpdf_nd->IsA() == TH2D::Class()) {
        hpdf[0] = dynamic_cast<TH2D*>(hpdf_nd)->ProjectionX("hpdf_x");
        hpdf[1] = dynamic_cast<TH2D*>(hpdf_nd)->ProjectionY("hpdf_y");
      }
      else if (hpdf_nd->IsA() == TH3D::Class()) {
        hpdf[0] = dynamic_cast<TH3D*>(hpdf_nd)->ProjectionX("hpdf_x");
        hpdf[1] = dynamic_cast<TH3D*>(hpdf_nd)->ProjectionY("hpdf_y");
        hpdf[2] = dynamic_cast<TH3D*>(hpdf_nd)->ProjectionZ("hpdf_y");
      }

      std::string n = signals[j].name;
      for (size_t k=0; k<observables.size(); k++) {
        hpdf[k]->SetLineColor(color[j]);
        if (fit_total[k] == NULL) {
          std::string hfname = "fit_total_" + signals[j].name;
          fit_total[k] = (TH1D*) hpdf[k]->Clone(hfname.c_str());
        }
        else {
          if (hpdf[k] && hpdf[k]->Integral() > 0) {
            fit_total[k]->Add(hpdf[k]);
          }
        }
  	TCanvas c1;
	fit_total[k]->Draw();
	c1.SaveAs((signals[j].name + ".pdf").c_str());


        if (n == "av_tl208" || n == "av_bi214" ||
            n == "water_tl208" || n == "water_bi214" ||
            n == "int_ropes_tl208" || n == "int_ropes_bi214" ||
            n == "hd_ropes_tl208" || n == "hd_ropes_bi214" ||
            n == "pmt_bg") {
          plots_external[k].add(hpdf[k], signals[j].title, "hist");
          std::string hname = "et" + signals[j].name + observables[k].name;
          if (external_total[k] == NULL) {
            external_total[k] = (TH1D*) hpdf[k]->Clone(hname.c_str());
          }
          else {
            if (hpdf[k] && hpdf[k]->Integral() > 0) {
              external_total[k]->Add((TH1D*) hpdf[k]->Clone(hname.c_str()));
            }
          }
        }
        else if (n != "zeronu" && n != "twonu" && n != "int_tl208" &&
                 n != "int_bi214" && n != "b8") {
          plots_cosmogenic[k].add(hpdf[k], signals[j].title, "hist");
          std::string hname = "ct" + signals[j].name + observables[k].name;
          if (cosmogenic_total[k] == NULL) {
            cosmogenic_total[k] = (TH1D*) hpdf[k]->Clone(hname.c_str());
          }
          else {
            cosmogenic_total[k]->Add((TH1D*) hpdf[k]->Clone(hname.c_str()));
          }
        }
        else {
          if (hpdf[k] && hpdf[k]->Integral() > 0) {
            plots_full[k].add(hpdf[k], signals[j].title, "hist");
          }
        }
      }
    }

    TCanvas c2;
    for (size_t j=0; j<observables.size(); j++) {
      TH1D* hdata = \
        (TH1D*) SpectralPlot::make_like(plots_full[j].histograms[0], "hdata");
      hdata->SetMarkerStyle(20);
      hdata->SetLineColor(kBlack);

      for (size_t idata=0; idata<data.size() / observables.size(); idata++) {
        hdata->Fill(data[idata * observables.size() + j]);
      }

      if (external_total[j] != NULL) {
        external_total[j]->SetLineColor(kOrange + 1);
        TH1D* et = (TH1D*) external_total[j]->Clone("et");
        et->SetLineStyle(2);
        plots_external[j].add(et, "Total", "hist");
        plots_full[j].add(external_total[j], "External", "hist");
      }
      if (cosmogenic_total[j] != NULL) {
        cosmogenic_total[j]->SetLineColor(kAzure + 1);
        TH1D* ct = (TH1D*) cosmogenic_total[j]->Clone("ct");
        ct->SetLineStyle(2);
        plots_cosmogenic[j].add(ct, "Total", "hist");
        plots_full[j].add(cosmogenic_total[j], "Cosmogenic", "hist");
      }
      if (fit_total[j] != NULL) {
        fit_total[j]->SetLineColor(kRed);
        plots_full[j].add(fit_total[j], "Fit", "hist");
      }

      plots_full[j].add(hdata, "Fake Data");

      plots_full[j].save(observables[j].name + "_spectrum_full.pdf");
      plots_external[j].save(observables[j].name + "_spectrum_external.pdf");
      plots_cosmogenic[j].save(observables[j].name +
                               "_spectrum_cosmogenic.pdf");
    }

    delete lspace;
    delete[] params_branch;
    delete[] params_fit;
  }

  hcontcoverage.Divide(&hcounts);
  hfccoverage.Divide(&hcounts);

  TCanvas ccov;
  hcontcoverage.SetLineColor(kRed);
  hcontcoverage.SetLineWidth(2);
  hcontcoverage.Draw();
  hcontcoverage.GetXaxis()->SetRangeUser(0, 20);
  hfccoverage.SetLineColor(kBlue);
  hfccoverage.SetLineWidth(2);
  hfccoverage.SetLineStyle(2);
  hfccoverage.Draw("same");
  ccov.SaveAs("coverage.pdf");
  ccov.SaveAs("coverage.C");

  return limits;
}


/**
 * Run the sensitivity calculation
 *
 * \param argc Argument count: must be 2
 * \param argv Arguments: executable path and configuration filename
 * \returns 0 on success, non-zero on failure
 */
int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " fit_configuration.json" << std::endl;
    exit(1);
  }

  // replace gRandom with something a little faster
  delete gRandom;
  gRandom = new TRandom2(0);
  gRandom->SetSeed(0);

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);

  // load configuration from json file
  std::string config_filename = std::string(argv[1]);
  FitConfig fc(config_filename);
  fc.print();

  // run ensemble
  std::map<std::string, TH1F> limits = \
    ensemble(fc.signals, fc.systematics, fc.observables, fc.cuts, fc.steps,
             fc.burnin_fraction, fc.signal_name, fc.signal_eff, fc.confidence,
             fc.experiments, fc.live_time, fc.debug_mode);

  // plot distribution of limits
  TCanvas c1;
  std::map<std::string, TH1F>::iterator it;
  for (it=limits.begin(); it!=limits.end(); it++) {
    it->second.DrawNormalized();
    c1.SaveAs((fc.output_file + "_limits_" + it->first + ".pdf").c_str());
    c1.SaveAs((fc.output_file + "_limits_" + it->first + ".C").c_str());

    // find median
    double xq[1] = {0.5};
    double yq[1];
    it->second.GetQuantiles(1, yq, xq);

    float lower_boundary, upper_boundary;
    bool is_limit;
    find_interval(&it->second, yq[0], lower_boundary, upper_boundary,
                  is_limit, false, 0.682);
    float lower_error = yq[0] - lower_boundary;
    float upper_error = upper_boundary - yq[0];

    std::cout << "Average-limit sensitivity " << it->first << " at "
              << fc.confidence * 100 << "\%: " << yq[0]
              << " -" << lower_error
              << " +" << upper_error
              << std::endl;
  }

  return 0;
}

