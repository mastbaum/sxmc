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
#include <utility>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
#include <TStyle.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TRandom.h>
#include <TRandom2.h>
#include <TMath.h>
#include <TError.h>

#include <sxmc/config.h>
#include <sxmc/generator.h>
#include <sxmc/mcmc.h>
#include <sxmc/utils.h>
#include <sxmc/likelihood.h>
#include <sxmc/plots.h>

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
 * \param confidence Desired confidence level for limits
 * \param nexperiments Number of fake experiments to run
 * \param live_time Experiment live time in years
 * \param debug_mode If true, accept and save all steps
 * \returns A list of the upper limits
 */
std::vector<float> ensemble(std::vector<Signal>& signals,
                             std::vector<Systematic>& systematics,
                             std::vector<Observable>& observables,
                             std::vector<Observable>& cuts,
                             unsigned steps, float burnin_fraction,
                             float confidence, unsigned nexperiments,
                             float live_time, const bool debug_mode,
                             std::string output_path) {
  std::vector<float> limits;

  for (size_t i=0;i<signals.size();i++){
    TFile f1((output_path+signals[i].name+"_pdf.root").c_str(),"RECREATE");
    TH1* hist = dynamic_cast<pdfz::EvalHist*> (signals[i].histogram)->DefaultHistogram();
    hist->Write();
    f1.Close();
  }

  for (unsigned i=0; i<nexperiments; i++) {
    std::cout << "Experiment " << i + 1 << " / " << nexperiments << std::endl;

    std::vector<double> params;
    std::vector<std::string> param_names;
    for (size_t j=0; j<signals.size(); j++) {
      params.push_back(signals[j].nexpected);
      param_names.push_back(signals[j].name);
    }
    for (size_t j=0; j<systematics.size(); j++) {
      params.push_back(systematics[j].mean);
      param_names.push_back(systematics[j].name);
    }


    // Make fake data
    std::pair<std::vector<float>, std::vector<int> > data = \
      make_fake_dataset(signals, systematics, observables, params, true);

    // Run MCMC
    MCMC mcmc(signals, systematics, observables);
    LikelihoodSpace* ls = mcmc(data.first, data.second, steps, burnin_fraction, debug_mode);

    // Write out samples for debugging
    TFile f((output_path+"lspace.root").c_str(), "recreate");
    TNtuple* lsclone = (TNtuple*) ls->GetSamples()->Clone("ls");
    lsclone->Write();
    lsclone->Delete();
    f.Close();

    ls->print_best_fit();
    ls->print_correlations();

    // Make spectral plots
    plot_fit(ls->get_best_fit(), live_time, signals,
             systematics, observables, data.first, data.second ,output_path);

    /*
    // Signal sensitivity, Bayesian for now
    TH1F* signal_projection = ls->get_projection(signal_name);
    int bin = 0;
    do {
      bin++;
    } while ((signal_projection->Integral(0, bin) /
              signal_projection->Integral()) < confidence);

    float limit = signal_projection->GetBinLowEdge(bin) +
                  signal_projection->GetBinWidth(bin);

    std::cout << "Signal limit: " << limit << std::endl;
    limits.push_back(limit);
    */

    delete ls;
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
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " fit_configuration.json output_path" << std::endl;
    exit(1);
  }

  // ROOT environment tweaks, including a faster RNG
  delete gRandom;
  gRandom = new TRandom2(0);
  gRandom->SetSeed(0);

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;

  std::string output_path = std::string(argv[2]);
  if (output_path.at( output_path.length() - 1 ) != '/'){
    output_path = output_path + "/";
  }

  // Load configuration from JSON file
  std::string config_filename = std::string(argv[1]);
  FitConfig fc(config_filename);
  fc.print();

  // Run ensemble
  std::vector<float> limits = \
    ensemble(fc.signals, fc.systematics, fc.observables, fc.cuts, fc.steps,
             fc.burnin_fraction, fc.confidence,
             fc.experiments, fc.live_time, fc.debug_mode, output_path);


  /*
  // Limits
  float average_limit;
  std::sort(limits.begin(), limits.end());
  float half = 1.0 * ((limits.size() + 1) / 2) - 1;
  average_limit = (limits.at(static_cast<int>(std::floor(half))) +
                   limits.at(static_cast<int>(std::ceil(half)))) / 2;

  std::cout << "Average-limit sensitivity " << average_limit << " at "
            << fc.confidence * 100 << "\% CL" << std::endl;

  TCanvas c1;
  TNtuple ntlimits("ntlimits", "ntlimits", "limit");
  for (size_t i=0; i<limits.size(); i++) {
    ntlimits.Fill(limits[i]);
  }
  ntlimits.Draw("limit>>_hlim", "", "goff");
  TH1F* hlim = dynamic_cast<TH1F*>(gDirectory->FindObject("_hlim"));
  assert(hlim);
  hlim->Draw();
  c1.SaveAs((fc.output_file + "_limits.pdf").c_str());
  c1.SaveAs((fc.output_file + "_limits.C").c_str());
  */

  return 0;
}

