/**
 * \file sensitivity.cpp
 * \brief Calculate sensitivity with an MCMC
 *
 * An ensemble of fake experiments is generated and each fit with an MCMC, with
 * parameters defined in the JSON configuration file. Limits are calculated for
 * each.
 */

#include <iostream>
#include <string>
#include <algorithm>
#include <assert.h>
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
 * \param e_range Energy range to consider
 * \param r_range Radial range to consider
 * \param steps Number of MCMC random walk steps to take
 * \param burnin_fraction Fraction of initial MCMC steps to throw out
 * \param signal_name The name of the Signal that is the signal
 * \param confidence Desired confidence level for limits
 * \param nexperiments Number of fake experiments to run
 * \returns A histogram of limits from the experiments
 */
TH1F* ensemble(std::vector<Signal>& signals, Range<float>& e_range,
               Range<float>& r_range, unsigned steps,
               float burnin_fraction, std::string signal_name,
               float confidence, unsigned nexperiments) {
  TH1F* limits = new TH1F("limits", "Background-fluctuation sensitivity",
                          300, 0, 300);

  for (unsigned i=0; i<nexperiments; i++) {
    std::cout << "Experiment " << i << " / " << nexperiments << std::endl;
    // make fake data
    FakeDataGenerator gen(signals, e_range, r_range);
    float* norms = new float[signals.size()];
    for (size_t i=0; i<signals.size(); i++) {
      norms[i] = signals[i].nexpected;
    }
    TCanvas c2;
    c2.SetLogz();
    TNtuple* data = gen(norms);
    data->Draw("e:r","","col z");
    c2.SaveAs("ee.pdf");

    // run mcmc
    MCMC mcmc(signals, data);
    TNtuple* lspace = mcmc(steps, burnin_fraction);

    // calculate signal sensitivity
    TH1F hproj("hproj", "hproj", 1000, 0, 100);
    lspace->Draw((signal_name + ">>hproj").c_str(), "", "goff");
    double limit = upper_limit(&hproj, confidence);
    std::cout << confidence * 100 << "\% limit: " << limit << std::endl;
    limits->Fill(limit);

    TCanvas c1;
    hproj.Rebin(10);
    hproj.Draw();
    c1.SaveAs("lproj.pdf");
    TFile f("lspace.root", "recreate");
    lspace->Write();
    f.Close();

    delete data;
    delete lspace;
    delete[] norms;
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

  // load configuration from json file
  std::string config_filename = std::string(argv[1]);
  FitConfig fc(config_filename);
  fc.print();

  // run ensemble
  TH1F* limits = ensemble(fc.signals, fc.e_range, fc.r_range, fc.steps,
                          fc.burnin_fraction, fc.signal_name,
                          fc.confidence, fc.experiments);

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

