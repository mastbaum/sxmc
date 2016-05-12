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

/**
 * \struct Observable
 *
 * A container for observable metadata
*/
struct Observable {
  std::string name;  //!< Name of the observable
  std::string title;  //!< Title in ROOT LaTeX format, for display
  std::string field;  //!< Name of the field (e.g. "energy")
  std::string units;  //!< Units as string, used in plotting
  bool logscale;  //!< Use log scale y in plots
  std::vector<float> yrange;  //!< y axis range for plots
  size_t field_index;  //!< Index in the sampled data for this field
  size_t bins;  //!< Number of bins
  float lower;  //!< Lower physical bound
  float upper;  //!< Upper physical bound
  float exclude_min;  //!< Minimum of excluded window
  float exclude_max;  //!< Maximum of excluded window
  bool exclude;  //!< Exclude a window inside the fit range
};


/**
 * \struct Systematic
 *
 * A container for systematic parameter metadata
*/
struct Systematic {
  std::string name;  //!< Name of the systematic parameter
  std::string title;  //!< Title in ROOT LaTeX format, for display
  std::string observable_field;  //!< Name of the field of the observable
  std::string truth_field;  //! Name of the field of the truth value
  size_t observable_field_index;  //!< Index of the observable field in the data
  size_t truth_field_index;  //!< Index of the truth field in the data
  pdfz::Systematic::Type type;  //!< The type of systematic
  size_t npars;  //!< Number of parameters in power series
  double* means;  //!< Mean values (power series)
  double* sigmas;  //!< Standard deviations
  bool fixed;  //! Fix the value of the parameter to the mean
};


/**
 * \class Signal
 *
 * A container for signal metadata and PDFs
*/
class Signal {
public:
  /**
   * Construct a Signal from a list of ROOT files.
   *
   * \param _name - A string identifier
   * \param _title - A title for plotting (ROOT LaTeX)
   * \oaram _nexpected - The expected number of events
   * \param _sigma - A Gaussian constraint on nexpected
   * \param _category - A category name, used to group signals for plotting
   * \param sample_fields - The names of the fields for the data samples
   * \param observables - A list of observables used in the fit
   * \param cuts - A list of observables to be applied as cuts
   * \param systematics - A list of systematics to be used in the fit
   * \param filenames - A list of ROOT file names with the PDF data
   * \param fixed - Normalization is fixed
  */
  Signal(std::string _name, std::string _title, float _nexpected,
         float _sigma, std::string _category,
         std::vector<std::string>& sample_fields,
         std::vector<Observable>& observables,
         std::vector<Observable>& cuts,
         std::vector<Systematic>& systematics,
         std::vector<std::string>& filenames,
         bool fixed=false);

  /**
   * Construct a Signal from a list of samples and weights.
   *
   * \param _name - A string identifier
   * \param _title - A title for plotting (ROOT LaTeX)
   * \oaram _nexpected - The expected number of events
   * \param _sigma - A Gaussian constraint on nexpected
   * \param _category - A category name, used to group signals for plotting
   * \param observables - A list of observables used in the fit
   * \param cuts - A list of observables to be applied as cuts
   * \param systematics - A list of systematics to be used in the fit
   * \param samples - The vector of data samples for the PDF
   * \param sample_fields - The names of the fields for the data samples
   * \param weights - The vector of sample weights for the PDF
   * \param fixed - Normalization is fixed
   */ 
  Signal(std::string _name, std::string _title, float _nexpected,
         float _sigma, std::string _category,
         std::vector<Observable>& observables,
         std::vector<Observable>& cuts,
         std::vector<Systematic>& systematics,
         std::vector<float>& samples,
         std::vector<std::string>& sample_fields,
         std::vector<int>& weights,
         bool fixed=false);

  /**
   * Get samples and weights as a pair of arrays.
  */
  //std::pair<std::vector<float>, std::vector<int> > get_samples();

  std::string name;  //!< String identifier
  std::string title;  //!< Histogram title in ROOT-LaTeX format
  std::string category;  //!< Category like external, cosmogenic, etc.
                         //!< for plotting purposes
  double nexpected;  //!< Events expected in this fit
  double sigma;  //!< Fractional uncertainty
  double efficiency;  //!< Fraction of generated events that make it past
                      //!< cuts (not counting the efficiency correction)
  bool fixed;  //!< Fix this parameter in the fit
  size_t nevents;  //!< Number of events in PDF
  size_t n_mc;  //!< Number of simulated events used to make pdf
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
  void build_pdfz(std::vector<float>& samples, std::vector<int>& weights,
                  int nfields, std::vector<Observable>& observables,
                  std::vector<Systematic>& systematics);

  /**
   * Compute an efficiency based on central value of systematics.
   *
   * \param systematics - A list of systematics to be used in the fit
   */
  void set_efficiency(std::vector<Systematic>& systematics);

public:
  /**
   * Apply an exclusion region that removes part of the dataset (with weights).
   *
   * Used e.g. to remove a blinded signal region.
   *
   * \param samples - The vector of data samples for the PDF
   * \param sample_fields - The names of the fields for the data samples
   * \param weights - The vector of sample weights for the PDF
   * \param observables - A list of observables used in the fit
   */
  static void apply_exclusions(std::vector<float>& samples,
                               std::vector<std::string>& sample_fields,
                               std::vector<int>& weights,
                               std::vector<Observable>& observables);

  /**
   * Apply an exclusion region that removes part of the dataset.
   *
   * Used e.g. to remove a blinded signal region.
   *
   * \param samples - The vector of data samples for the PDF
   * \param sample_fields - The names of the fields for the data samples
   * \param observables - A list of observables used in the fit
   */
  static void apply_exclusions(std::vector<float>& samples,
                               std::vector<std::string>& sample_fields,
                               std::vector<Observable>& observables) {
    std::vector<int> fake;
    apply_exclusions(samples, sample_fields, fake, observables);
  };

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
                                      std::vector<std::string>& sample_fields,
                                      std::vector<std::string>& dataset_fields,
                                      std::vector<Observable>& cuts);
};

#endif  // __SIGNALS_H__

