#ifndef __CONFIG_H__
#define __CONFIG_H__

/**
 * \file config.h
 *
 * Fit configuration management and config file parsing.
*/

#include <string>
#include <vector>
#include <json/value.h>

#include <sxmc/utils.h>
#include <sxmc/signals.h>

/**
 * \struct SignalParams
 *
 * Container for signal parameters extracted from config JSON object.
*/
class SignalParams {
public:
  /** Default constructor, does no initialization. */
  SignalParams() {}

  /**
   * Extract the parameters from a JSON::Value configuration object.
   *
   * \param params The JSON object
   * \param scale Scaling factor applied to the expected rate and sigma
  */
  SignalParams(const Json::Value& params, float scale=1.0);

  float nexpected;  //!< Number of events expected
  float sigma;  //!< Gaussian constraint
  std::string title;  //!< Title of signal (for plotting)
  std::string category;  //!< Category (for plotting)
  std::vector<std::string> files;  //!< List of filenames with dataset
};


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
  virtual ~FitConfig() {
    delete data;
  }

  /** Pretty-print the fit parameters */
  void print() const;

  unsigned experiments;  //!< Number of experiments in ensemble
  unsigned steps;  //!< Number of MCMC steps
  float confidence;  //!< Confidence level for results (e.g. 0.9)
  float live_time;  //!< Experiment live time in years
  float efficiency_correction;  //!< Overall efficiency correction
  float burnin_fraction;  //!< Fraction of steps to use for burn-in period
  bool debug_mode;  //!< Enable/disable debugging mode (accept/save all)
  std::string output_file;  //!< Base filename for output
  std::string signal_name;  //!< Name of the signal of interest (if any)
  std::vector<Signal> signals;  //!< Signal histograms and metadata
  std::vector<Systematic> systematics;  //!< Systematics used in PDFs
  std::vector<Observable> observables;  //!< Observables used in PDFs
  std::vector<Observable> cuts;  //!< Cuts applied before fit
  std::vector<std::string> sample_fields;  //!< Names of sample array fields:
                                           //!< observables and systematics
  Signal* data;  //!< Data to fit, if any
};

#endif  // __CONFIG_H__

