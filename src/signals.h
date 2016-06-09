#ifndef __SIGNALS_H__
#define __SIGNALS_H__

/**
 * \file signals.h
 *
 * Signal-related data structures.
 *
 * (Not to be confused with signal.h, which defines UNIX signals.)
 */

#include <string>
#include <vector>
#include <sxmc/pdfz.h>

namespace Json {
  class Value;
}

/**
 * \struct Observable
 *
 * A container for observable metadata
*/
struct Observable {
  Observable() {}

  // Load from JSON configuration
  Observable(const std::string _name, const Json::Value& config);

  std::string name;  //!< Name of the observable
  std::string title;  //!< Title in ROOT LaTeX format, for display
  std::string field;  //!< Name of the field (e.g. "energy")
  std::string units;  //!< Units as string, used in plotting
  std::vector<float> yrange;  //!< y axis range for plots
  size_t field_index;  //!< Index in the sampled data for this field
  size_t bins;  //!< Number of bins
  float lower;  //!< Lower physical bound
  float upper;  //!< Upper physical bound
  bool logscale;  //!< Use log scale y in plots
};


/**
 * \struct Systematic
 *
 * A container for systematic parameter metadata
*/
struct Systematic {
  Systematic() {}

  // Load from JSON configuration
  Systematic(const std::string _name, const Json::Value& config);

  std::string name;  //!< Name of the systematic parameter
  std::string title;  //!< Title in ROOT LaTeX format, for display
  std::string observable_field;  //!< Name of the field of the observable
  std::string truth_field;  //! Name of the field of the truth value
  double* means;  //!< Mean values (power series)
  double* sigmas;  //!< Standard deviations
  size_t observable_field_index;  //!< Index of the observable field in the data
  size_t truth_field_index;  //!< Index of the truth field in the data
  size_t npars;  //!< Number of parameters in power series
  std::vector<short> pidx;  //!< Global index for pdfz parameter array offsetting
  pdfz::Systematic::Type type;  //!< The type of systematic
  bool fixed;  //! Fix the value of the parameter to the mean
};


///**
// * \struct Rate
// *
// * A container for a signal rate.
//*/
//struct Rate {
//  Rate() {}
//
//  // Load from JSON configuration
//  Rate(const std::string& _name, const Json::Value& params);
//
//  std::string name;  //!< String identifier
//  size_t index;  //!< Index in the list of rates
//  float rate;  //!< Number of events expected
//  float sigma;  //!< Gaussian constraint
//  bool fixed;  //!< Rate is fixed
//};


/**
 * \class Signal
 *
 * A container for signal metadata and PDFs
*/
class Signal {
public:
  Signal() {}

  /**
   * Construct a Signal from a list of ROOT files.
   *
   * \param _name - A string identifier
   * \param _title - A title for plotting (ROOT LaTeX)
   * \oaram _nexpected - The expected number of events
   * \param _sigma - A Gaussian constraint on nexpected
   * \param sample_fields - The names of the fields for the data samples
   * \param observables - A list of observables used in the fit
   * \param cuts - A list of observables to be applied as cuts
   * \param systematics - A list of systematics to be used in the fit
   * \param filenames - A list of ROOT file names with the PDF data
   * \param fixed - Normalization is fixed
  */
  Signal(std::string _name, std::string _title,
         std::string _filename, unsigned _dataset,
         double _nexpected, double _sigma, bool _fixed,
         std::vector<std::string>& sample_fields,
         std::vector<Observable>& observables,
         std::vector<Observable>& cuts,
         std::vector<Systematic>& systematics);

  std::string name;  //!< String identifier
  std::string title;  //!< Histogram title in ROOT-LaTeX format
  std::string filename;
  unsigned dataset;  //!< Dataset identifier (cf. multi-phase fitting)
  double nexpected;  //!< Events expected in this fit
  float sigma;  //!< Gaussian constraint
  bool fixed;  //!< Rate is fixed
  size_t n_mc;  //!< Number of simulated events used to make pdf
  double nevents;  //!< Events in the PDF
  double efficiency;  //!< Fraction of generated events that make it past
                      //!< cuts (not counting the efficiency correction)
  pdfz::Eval* histogram;  //!< PDF

protected:
  /**
   * Construct the pdfz histogram object.
   *
   * \param samples - The vector of data samples for the PDF
   * \param weights - The vector of sample weights for the PDF
   * \param nfields - The number of fields (columns) in the sample array
   * \param observables - A list of observables used in the fit
   * \param systematics - A list of systematics to be used in the fit
   */
  void build_pdfz(std::vector<float>& samples, int nfields,
                  std::vector<Observable>& observables,
                  std::vector<Systematic>& systematics);

  /**
   * Compute an efficiency based on central value of systematics.
   *
   * \param systematics - A list of systematics to be used in the fit
   */
  void set_efficiency(std::vector<Systematic>& systematics);

public:
  /**
   * Copy the requested sample fields from the fields of the same name
   * in the dataset array.
   *
   * \param samples - The vector of data samples for the PDF (by reference)
   * \param dataset - The dataset to load in, as a float array
   * \param sample_fields - The names of the fields for the samples array
   * \param sample_fields - The names of the fields in the dataset array
   * \param cuts - A list of observables to be applied as cuts
   */
  static void read_dataset_to_samples(std::vector<float>& samples,
                                      std::vector<float>& dataset,
                                      unsigned dataset_id,
                                      std::vector<std::string>& sample_fields,
                                      std::vector<std::string>& dataset_fields,
                                      std::vector<Observable>& cuts);
};

#endif  // __SIGNALS_H__

