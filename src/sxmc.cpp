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
 * \param signal_name The name of the Signal that is the signal
 * \param signal_eff Efficiency for detecting signal
 * \param confidence Desired confidence level for limits
 * \param nexperiments Number of fake experiments to run
 * \param live_time Experiment live time in years
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
                                     float live_time, const bool debug_mode,
                                     std::string output_path)
{
  for (size_t i=0;i<signals.size();i++){
    TFile f1((output_path+signals[i].name+"_pdf.root").c_str(),"RECREATE");
    TH1* hist = dynamic_cast<pdfz::EvalHist*> (signals[i].histogram)->DefaultHistogram();
    hist->Write();
    f1.Close();
  }

  std::map<std::string, TH1F> limits;
  //limits["counts_proj"] = TH1F("counts_proj", ";Counts in fit range;Fraction",
  //                             10000, 0, 500);
  //limits["counts_cont"] = TH1F("counts_cont", ";Counts in fit range;Fraction",
  //                             10000, 0, 500);

  for (unsigned i=0; i<nexperiments; i++) {
    std::cout << "Experiment " << i + 1 << " / " << nexperiments << std::endl;

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


    // Make fake data
    std::vector<float> data = \
      make_fake_dataset(signals, systematics, observables, params, true);

    // Run MCMC
    MCMC mcmc(signals, systematics, observables);
    LikelihoodSpace* ls = mcmc(data, steps, burnin_fraction, debug_mode);

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
             systematics, observables, data,output_path);

    // Signal sensitivity

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
  std::map<std::string, TH1F> limits = \
    ensemble(fc.signals, fc.systematics, fc.observables, fc.cuts, fc.steps,
             fc.burnin_fraction, fc.signal_name, fc.signal_eff, fc.confidence,
             fc.experiments, fc.live_time, fc.debug_mode, output_path);


  // Plot distribution of limits
/*
  TCanvas c1;
  std::map<std::string, TH1F>::iterator it;
  for (it=limits.begin(); it!=limits.end(); it++) {
    it->second.DrawNormalized();
    c1.SaveAs((fc.output_file + "_limits_" + it->first + ".pdf").c_str());
    c1.SaveAs((fc.output_file + "_limits_" + it->first + ".C").c_str());
*/

/*
    // Find median
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
*/

  return 0;
}

