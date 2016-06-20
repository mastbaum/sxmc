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
#include <TStyle.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TRandom.h>
#include <TRandom2.h>
#include <TROOT.h>
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
 * Run experiments, fit each with MCMC.
 *
 * \param fc A FitConfig defining the fit parameters
 * \param output_path Location to output files
 * \returns A list of the upper limits
 */
std::vector<float> ensemble(FitConfig& fc, std::string output_path) {
  if (fc.plots && fc.samples == "") {
    for (size_t i=0; i<fc.signals.size(); i++) {
      const char* output_file = \
        (output_path + fc.signals[i].name + "_pdf.root").c_str();
      TFile f1(output_file, "RECREATE");
      TH1* hist = \
        dynamic_cast<pdfz::EvalHist*>(fc.signals[i].histogram)->DefaultHistogram();
      hist->Write();
      f1.Close();
    }
  }

  std::vector<float> limits;

  for (unsigned i=0; i<fc.nexperiments; i++) {
    std::cout << "Experiment " << i + 1
              << " / " << fc.nexperiments << std::endl;

    // Make fake data or load a dataset
    std::vector<float> samples;
    if (fc.data.empty()) {
      std::cout << "ensemble: Sampling fake dataset " << i << std::endl;
      samples = \
        make_fake_dataset(fc.signals, fc.systematics, fc.observables, true);
    }
    else {
      for (std::map<unsigned, std::vector<Signal> >::iterator it=fc.data.begin();
           it!=fc.data.end(); ++it) {
        std::cout << "ensemble: Loading dataset "
                  << it->first << "." << i << " ("
                  << it->second[i].filename << ")" << std::endl;

        dynamic_cast<pdfz::EvalHist*>(it->second[i].histogram)->GetSamples(samples);
      }
    }

    // Run MCMC (or load MCMC samples from a file)
    TFile* sample_file = NULL;
    LikelihoodSpace* ls = NULL;
    if (fc.samples != "") {
      std::cout << "ensemble: Loading samples from "
                << fc.samples << std::endl;
      sample_file = TFile::Open(fc.samples.c_str(),"update");
      assert(sample_file && sample_file->IsOpen());
      TNtuple* nt = (TNtuple*) sample_file->Get("ls")->Clone("_ls");
      assert(nt && nt->GetEntries() > 0);
      std::cout << "ensemble: Loaded " << nt->GetEntries() << " samples"
                << std::endl;
      ls = new LikelihoodSpace(nt, fc.confidence, fc.error_type);
    }
    else {
      MCMC mcmc(fc.sources, fc.signals, fc.systematics, fc.observables);
      ls = mcmc(samples, fc.nsteps, fc.burnin_fraction, fc.debug_mode);
    }

    ls->print_best_fit();
    ls->print_correlations();

    // Make projection plots
    if (fc.plots) {
      plot_fit(ls->get_best_fit(), 1.0, fc.sources, fc.signals, fc.systematics,
               fc.observables, fc.datasets, samples, output_path);
    }

    // Build a list of upper limits
    std::map<std::string, Interval> best_fit = ls->get_best_fit();

    if (fc.signal_name != "" &&
        best_fit.find(fc.signal_name) != best_fit.end()) {
      Interval bfi = best_fit[fc.signal_name];
      std::cout << "ensemble: Signal "
                << fc.signal_name << ": " << bfi.str() << std::endl;

      if (!bfi.one_sided) {
        std::cerr << "ensemble: Warning: Two-sided limit!" << std::endl;
      }

      std::cout << "ensemble: lower = " << bfi.lower << ", "
                << "upper = " << bfi.upper << ", "
                << "coverage = " << bfi.coverage << std::endl;

      limits.push_back(bfi.upper);
    }

    // Write out samples
    if (!sample_file) {
      std::ostringstream lsfile;
      lsfile << output_path << fc.output_prefix << "_" << i << ".root";
      TFile f(lsfile.str().c_str(), "recreate");
      assert(f.IsOpen());
      TNtuple* lsclone = \
        dynamic_cast<TNtuple*>(ls->get_samples()->Clone("ls"));
      assert(lsclone);
      lsclone->Write();
      lsclone->Delete();
      f.Close();
    }

    delete sample_file;
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
  gErrorIgnoreLevel = kError;
  gROOT->SetBatch(true);

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
  std::vector<float> limits = ensemble(fc, output_path);

  std::cout << "sxmc: Upper limits: ";
  for (size_t i=0; i<limits.size(); i++) {
    std::cout << limits[i] << ", ";
  }
  std::cout << std::endl;
  std::cout << "sxmc: Median upper limit: "
            << (!limits.empty() ? median(limits) : -1) << std::endl;

  return 0;
}

