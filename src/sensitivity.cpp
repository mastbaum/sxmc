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
#include <algorithm>
#include <assert.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TRandom.h>
#include <TRandom2.h>
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
    integral = h->Integral(thisbin--, h->GetNbinsX());
  }
  return h->GetBinLowEdge(thisbin);
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
 * \param steps Number of MCMC random walk steps to take
 * \param burnin_fraction Fraction of initial MCMC steps to throw out
 * \param signal_name The name of the Signal that is the signal
 * \param confidence Desired confidence level for limits
 * \param nexperiments Number of fake experiments to run
 * \returns A histogram of limits from the experiments
 */
TH1F* ensemble(std::vector<Signal>& signals,
               std::vector<Systematic>& systematics,
               std::vector<Observable>& observables,
               unsigned steps, float burnin_fraction, std::string signal_name,
               float confidence, unsigned nexperiments) {
  TH1F* limits = new TH1F("limits", "Background-fluctuation sensitivity",
                          100, 0, 100);

  for (unsigned i=0; i<nexperiments; i++) {
    std::cout << "Experiment " << i << " / " << nexperiments << std::endl;

    // make fake data
    std::vector<float> params;
    std::vector<std::string> param_names;
    for (size_t j=0; j<signals.size(); j++) {
      params.push_back(signals[j].nexpected);
      param_names.push_back(signals[j].name);
    }
    for (size_t j=0; j<systematics.size(); j++) {
      params.push_back(systematics[j].mean);
      param_names.push_back(systematics[j].name);
    }
    std::vector<float> data = make_fake_dataset(signals, systematics,
                                                observables, params, true);

    // run mcmc
    MCMC mcmc(signals, systematics, observables);
    TNtuple* lspace = mcmc(data, steps, burnin_fraction);

    // calculate signal sensitivity
    TH1F hproj("hproj", "hproj", 1000, 0, 100);
    lspace->Draw((signal_name + ">>hproj").c_str(), "", "goff");
    double limit = upper_limit(&hproj, confidence);
    std::cout << confidence * 100 << "\% limit: " << limit << std::endl;
    limits->Fill(limit);

    // extract likelihood-maximizing parameters
    float* params_branch = new float[params.size()];
    for (size_t j=0; j<params.size(); j++) {
      lspace->SetBranchAddress(param_names[j].c_str(), &params_branch[j]);
    }
    float ml_branch;
    lspace->SetBranchAddress("likelihood", &ml_branch);
    float* params_fit = new float[params.size()];
    float ml = 1e9;
    for (int j=0; j<lspace->GetEntries(); j++) {
      lspace->GetEntry(j);
      if (ml_branch < ml) {
        ml = ml_branch;
        for (size_t k=0; k<params.size(); k++) {
          params_fit[k] = params_branch[k];
        }
      }
    }

    std::cout << "-- Best fit: NLL = " << ml << " --" << std::endl;
    for (size_t j=0; j<params.size(); j++) {
      std::cout << " " << param_names[j] << ": " << params_fit[j]
                << " (" << params[j] <<  ")" << std::endl;
    }

    // plot and save this fit
    bool energy_observable_found = false;
    size_t energy_observable_index = 0;
    for (size_t j=0; j<observables.size(); j++) {
      if (observables[i].field == "e") {
        energy_observable_found = true;
        energy_observable_index = i;
        break;
      }
    }
    assert(energy_observable_found);

    SpectralPlot p_all;
    SpectralPlot p_external;
    SpectralPlot p_cosmo;

    pdfz::EvalHist* phist = \
      dynamic_cast<pdfz::EvalHist*>(signals[0].histogram);

    hemi::Array<unsigned> norms_buffer(signals.size(), true);
    norms_buffer.writeOnlyHostPtr();
    phist->SetNormalizationBuffer(&norms_buffer, 0);

    hemi::Array<float> param_buffer(params.size(), true);
    for (size_t j=0; j<params.size(); j++) {
      param_buffer.writeOnlyHostPtr()[j] = params_fit[j];
    }
    phist->SetParameterBuffer(&param_buffer, signals.size());

    TH1* hpdf0 = phist->CreateHistogram();
    TH1* hdata = SpectralPlot::make_like(hpdf0, "hdata");
    hdata->SetAxisRange(1e-1, 1e3, "Y");
    hdata->SetAxisRange(1.5, 5.0, "X");
    hdata->SetMarkerStyle(20);
    hdata->SetLineColor(kBlack);

    for (size_t idata=0; idata<data.size(); idata++) {
      hdata->Fill(data[idata * observables.size() + energy_observable_index]);
    }

    p_all.add(hdata, "Data");

    TH1* hsum = SpectralPlot::make_like(hpdf0, "hdata");
    hsum->SetLineColor(kRed);

    TH1* hcosmo = SpectralPlot::make_like(hpdf0, "hcosmo");
    hcosmo->SetLineColor(kAzure + 1);

    TH1* hexternal = SpectralPlot::make_like(hcosmo, "hexternal");
    hexternal->SetLineColor(kOrange + 1);

    for (size_t j=0; j<signals.size(); j++) {
      // create TH1 histogram from pdfz::EvalHist
      phist = dynamic_cast<pdfz::EvalHist*>(signals[j].histogram);
      phist->SetParameterBuffer(&param_buffer, signals.size());
      phist->SetNormalizationBuffer(&norms_buffer, j);
      TH1* hpdf = phist->CreateHistogram();

      TH1* hs = (TH1*) hpdf->Clone("hs");
      int bin1 = hs->FindBin(1.5);
      int bin2 = hs->FindBin(5.0);
      //float syst_rescale = 
      //1.0 * norms_buffer.readOnlyHostPtr()[i] / signals[i].nevents;
      float syst_rescale = 1.0;
      hs->Scale(params_fit[j] * syst_rescale / hpdf->Integral(bin1, bin2));
      hs->SetLineColor(color[j + 1]);
      hsum->Add(hs);

      hs->SetAxisRange(1e-1, 1e3, "Y");
      hs->SetAxisRange(1.5, 5.0, "X");

      std::string n = signals[j].name;
      if (n == "av_tl208" || n == "av_bi214" || n == "water_tl208" ||
          n == "water_bi214" || n == "pmt_bg") {
        p_external.add(hs, signals[j].title, "hist");
        hexternal->Add(hs);
      }
      else if (n != "zeronu" && n != "twonu" && n != "int_tl208" &&
               n != "int_bi214" && n != "b8") {
        p_cosmo.add(hs, signals[j].title, "hist");
        hcosmo->Add(hs);
      }
      else {
        p_all.add(hs, signals[j].title, "hist");
      }
      delete hs;
    }

    p_all.add(hexternal, "External", "hist");
    p_all.add(hcosmo, "Cosmogenics", "hist");
    p_all.add(hsum, "Fit", "hist");
    p_all.save("spectrum_all.pdf");

    hexternal->SetLineStyle(2);
    p_external.add(hexternal, "Sum", "hist");
    p_external.save("spectrum_external.pdf");

    hcosmo->SetLineStyle(2);
    p_cosmo.add(hcosmo, "Sum", "hist");
    p_cosmo.save("spectrum_cosmogenic.pdf");

    TCanvas c2;
    hproj.Rebin(10);
    hproj.Draw();
    c2.SaveAs("lproj.pdf");

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

    TFile f("lspace.root", "recreate");
    lspace->Write();
    f.Close();

    delete lspace;
    delete[] params_branch;
  }

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
  TH1F* limits = ensemble(fc.signals, fc.systematics, fc.observables, fc.steps,
                          fc.burnin_fraction, fc.signal_name, fc.confidence,
                          fc.experiments);

  // plot distribution of limits
  TCanvas c1;
  limits->Draw();
  c1.SaveAs((fc.output_file + "_limits.pdf").c_str());
  c1.SaveAs((fc.output_file + "_limits.C").c_str());

  // find median
  double xq[1] = {0.5};
  double yq[1];
  limits->GetQuantiles(1, yq, xq);

  std::cout << "Median " << fc.confidence * 100 << "\% limit: " << yq[0]
            << std::endl;

  delete limits;

  return 0;
}

