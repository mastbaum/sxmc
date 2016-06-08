/**
 * \file sxmc.cpp
 * \brief Calculate sensitivity with an MCMC
 *
 * An ensemble of fake experiments is generated and each fit with an MCMC, with
 * parameters defined in the JSON configuration file. Limits are calculated for
 * each.
 */

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <unistd.h>
#include <iostream>
#include <vector>
#include <utility>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TH1F.h>
#include <TRandom.h>
#include <TRandom2.h>
#include <TError.h>

#include <sxmc/config.h>
#include <sxmc/error_estimator.h>
#include <sxmc/generator.h>
#include <sxmc/mcmc.h>
#include <sxmc/utils.h>
#include <sxmc/likelihood.h>
#include <sxmc/plots.h>
#include <sxmc/utils.h>

/**
 * Run an ensemble of independent fake experiments
 *
 * Run experiments, fit each with MCMC.
 *
 * \param signals List of Signals defining PDFs, rates, etc.
 * \param systematics List of Systematics applied to PDFs
 * \param observables List of Observables common to PDFs
 * \param cuts List of cuts applied to data set
 * \param data Signal containing data to fit, or NULL to sample
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
                            std::vector<Signal>* data,
                            unsigned steps, float burnin_fraction,
                            float confidence, unsigned nexperiments,
                            float live_time, const bool debug_mode,
                            std::string output_path, std::string signal_name,
                            std::string prefix, bool plots) {
  std::vector<float> limits;

  if (plots) {
    for (size_t i=0; i<signals.size(); i++) {
      const char* output_file = \
        (output_path + signals[i].name + "_pdf.root").c_str();
      TFile f1(output_file, "RECREATE");
      TH1* hist = \
        dynamic_cast<pdfz::EvalHist*>(signals[i].histogram)->DefaultHistogram();
      hist->Write();
      f1.Close();
    }
  }

  for (unsigned i=0; i<nexperiments; i++) {
    std::cout << "Experiment " << i + 1 << " / " << nexperiments << std::endl;

    std::vector<double> params;
    std::vector<std::string> param_names;
    for (size_t j=0; j<signals.size(); j++) {
      params.push_back(signals[j].nexpected);
      param_names.push_back(signals[j].name);
    }

    std::vector<float> samples;
    std::vector<int> weights;
    if (!data) {
      // Make fake data
      std::cout << "ensemble: Sampling fake dataset " << i << std::endl;
      std::pair<std::vector<float>, std::vector<int> > fakedata = \
        make_fake_dataset(signals, systematics, observables, params, true);
      samples = fakedata.first;
      weights = fakedata.second;
    }
    else {
      std::cout << "ensemble: Using dataset " << i << std::endl;
      samples = \
        dynamic_cast<pdfz::EvalHist*>(data->at(i).histogram)->GetSamples();
      weights.resize(samples.size(), 1);
    }

    // Run MCMC
    MCMC mcmc(signals, systematics, observables);
    LikelihoodSpace* ls = \
      mcmc(samples, weights, steps, burnin_fraction, debug_mode);

    // Write out samples for debugging
    std::ostringstream lsfile;
    lsfile << output_path << prefix << "_" << i << ".root";
    TFile f(lsfile.str().c_str(), "recreate");
    TNtuple* lsclone = dynamic_cast<TNtuple*>(ls->get_samples()->Clone("ls"));
    lsclone->Write();
    lsclone->Delete();
    f.Close();

    ls->print_best_fit();
    ls->print_correlations();

    // Make projection plots
    if (plots) {
      plot_fit(ls->get_best_fit(), live_time, signals, systematics,
               observables, samples, weights, output_path);
    }

    // Build a list of upper limits
    float ml;
    std::map<std::string, Interval> best_fit = \
      ls->extract_best_fit(ml, confidence, ERROR_PROJECTION);

    if (signal_name != "" && best_fit.find(signal_name) != best_fit.end()) {
      Interval bfi = best_fit[signal_name];
      std::cout << "ensemble: Signal "
                << signal_name << ": " << bfi.str() << std::endl;

      if (!bfi.one_sided) {
        std::cerr << "ensemble: Warning: Two-sided limit!" << std::endl;
      }

      std::cout << "ensemble: lower = " << bfi.lower << ", "
                << "upper = " << bfi.upper << ", "
                << "coverage = " << bfi.coverage << std::endl;

      limits.push_back(bfi.upper);
    }

    delete ls;
  }

  return limits;
}


/**
 * Check the output path; create if necessary.
 *
 * \param path The path for the output directory
 * \returns 0 on success, non-zero on failure
 */
int check_create_output(std::string path) {
  struct stat st;
  if (stat(path.c_str(), &st) != 0) {
    if (mkdir(path.c_str(), 0777) != 0 && errno != EEXIST) {
      std::cerr << "Error creating output directory, errno = "
                << errno << std::endl;
      return -1;
    }
  }
  else if (!S_ISDIR(st.st_mode)) {
    std::cerr << "Output path exists and is not a directory" << std::endl;
    return -1;
  }
  return 0;
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
    std::cerr << "Usage: " << argv[0]
              << " fit_configuration.json output_path"
              << std::endl;
    exit(1);
  }

  // ROOT environment tweaks, including a faster RNG
  delete gRandom;
  gRandom = new TRandom2(0);
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;

  // Set and check/create output path
  std::string output_path = std::string(argv[2]);
  if (output_path.at(output_path.length() - 1) != '/'){
    output_path = output_path + "/";
  }
  if (check_create_output(output_path) != 0) {
    std::cerr << "Error with output path." << std::endl;
    throw(2);
  }

  // Load configuration from JSON file
  std::cout << "sxmc: Loading configuration..." << std::endl;
  std::string config_filename = std::string(argv[1]);
  std::cout << "sxmc: Configuration: " << config_filename << std::endl;
  FitConfig fc(config_filename);
  fc.print();

  gRandom->SetSeed(fc.seed);

  // Run ensemble
  std::cout << "sxmc: Running ensemble..." << std::endl;
  std::vector<float> limits = \
    ensemble(fc.signals, fc.systematics, fc.observables, fc.cuts, fc.data,
             fc.steps, fc.burnin_fraction, fc.confidence, fc.experiments,
             fc.live_time, fc.debug_mode, output_path, fc.signal_name,
             fc.prefix, fc.plots);

  std::cout << "sxmc: Upper limits: ";
  for (size_t i=0; i<limits.size(); i++) {
    std::cout << limits[i] << ", ";
  }
  std::cout << std::endl;
  std::cout << "sxmc: Median upper limit: "
            << (limits.size() > 0 ? median(limits) : -1) << std::endl;

  return 0;
}

