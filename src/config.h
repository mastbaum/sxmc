#ifndef __CONFIG_H__
#define __CONFIG_H__

/**
 * \file config.h
 *
 * Fit configuration management and config file parsing.
*/

#include <map>
#include <set>
#include <string>
#include <vector>

#include <sxmc/signals.h>

/**
 * \class FitConfig
 * \brief Manages the configuration of a fit
 *
 * A utility class for reading in configuration parameters from a JSON file
 * and storing the parameters of the fit.
*/
class FitConfig {
public:
  /**
   * Constructor
   *
   * \param filename Name of the JSON file to load the configuration from
  */
  FitConfig(std::string filename);

  /** Destructor */
  virtual ~FitConfig() {}

  /** Pretty-print the fit parameters */
  void print() const;

  unsigned nexperiments;  //!< Number of experiments in ensemble
  unsigned nsteps;  //!< Number of MCMC steps
  float confidence;  //!< Confidence level for results (e.g. 0.9)
  float burnin_fraction;  //!< Fraction of steps to use for burn-in period
  bool debug_mode;  //!< Enable/disable debugging mode (accept/save all)
  bool plots;  //!< Write out fit plots and PDF root files
  long seed;  //!< RNG seed
  std::string signal_name;  //!< Name of the signal of interest (if any)
  std::string output_prefix;  //!< Prefix for output ROOT file name
  std::vector<Signal> signals;  //!< Signal histograms and metadata
  std::vector<Systematic> systematics;  //!< Systematics used in PDFs
  std::vector<Observable> observables;  //!< Observables used in PDFs
  std::vector<Observable> cuts;  //!< Cuts applied before fit
  std::map<unsigned, std::vector<Signal> > data;  //!< Data to fit, if any
  std::vector<std::string> sample_fields;  //!< Names of sample array fields:
  std::set<unsigned> datasets;  //!< Dataset tags
};

#endif  // __CONFIG_H__

