#ifndef __CONFIG_H__
#define __CONFIG_H__

/**
 * \file config.h
 * \brief Fit configuration management and config file parsing.
*/

#include <map>
#include <set>
#include <string>
#include <vector>

#include <sxmc/error_estimator.h>
#include <sxmc/signal.h>
#include <sxmc/observable.h>
#include <sxmc/systematic.h>

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
  float burnin_fraction;  //!< Fraction of steps to use for burn-in periods
  float mc_scale;  //!< Global MC normalization scaling
  bool debug_mode;  //!< Enable/disable debugging mode (accept/save all steps)
  bool plots;  //!< Write out fit plots and PDF root files
  bool sensitivity;  //!< Sensitivity mode, no signal in fake data
  long seed;  //!< RNG seed
  std::string signal_name;  //!< Name of the signal of interest (if any)
  std::string output_prefix;  //!< Prefix for output ROOT file name
  std::vector<Signal> signals;  //!< Signal histograms and metadata
  std::vector<Systematic> systematics;  //!< Systematics used in PDFs
  std::vector<Observable> observables;  //!< Observables used in PDFs
  std::vector<Observable> cuts;  //!< Cuts applied before fit
  std::vector<Source> sources;  //!< Sources for correlated rates
  std::map<unsigned, std::vector<Signal> > data;  //!< Data to fit, if any
  std::vector<std::string> sample_fields;  //!< Names of sample array fields
  std::set<unsigned> datasets;  //!< Dataset tags
  std::string samples;  //!< Path to externally-supplied samples
  ErrorType error_type;  //!< Error calculation method
};

#endif  // __CONFIG_H__

