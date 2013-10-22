/**
 * \file sxmc.cpp
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
#include "likelihood.h"

/** A ROOT color palette in which sequential colors look okay. */
static const int color[28] = {kRed,      kGreen,    kBlue,      kMagenta,
                              kCyan,     kYellow,   kOrange,    kViolet+2,
                              kRed+2,    kGreen+2,  kBlue+2,    kMagenta+2,
                              kCyan+2,   kYellow+2, kOrange+2,  kRed-7,
                              kGreen-7,  kBlue-7,   kMagenta-7, kCyan-7,
                              kYellow-7, kOrange-7, kRed-6,     kAzure+1,
                              kTeal+1,   kSpring-9, kAzure-9};


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
    std::cout << "Experiment " << i + 1 << " / " << nexperiments << std::endl;

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
    TNtuple* ls_samples = mcmc(data, steps, burnin_fraction, debug_mode);

    TFile f("lspace.root", "recreate");
    TNtuple* lsclone = (TNtuple*) ls_samples->Clone("ls_samples");
    lsclone->Write();
    //lsclone->Delete();
    f.Close();

    // calculate signal sensitivity
    LikelihoodSpace lspace(ls_samples);
    lspace.print_best_fit();
    lspace.print_correlations();
/*
    // plot this fit
    std::vector<SpectralPlot> plots_full;
    std::vector<SpectralPlot> plots_external;
    std::vector<SpectralPlot> plots_cosmogenic;
    for (size_t j=0; j<observables.size(); j++) {
      Observable* o = &observables[j];
      std::stringstream ytitle;
      ytitle << "Counts/" << std::setprecision(3)
             << (o->upper - o->lower) / o->bins << " " << o->units
             << "/" << live_time << " y";
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
        // hack to include only radius ROI
        if (j == 0 && data[idata * observables.size() + 1] > radius_cut) {
          continue;
        }
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
*/
    //ls_samples->Delete();

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

/*
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
*/
  }

  return 0;
}

