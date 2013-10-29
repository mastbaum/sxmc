/**
 * \file signals.h
 *
 * Signal-related data structures.
 *
 * (Not to be confused with signal.h, which defines UNIX signals.)
 */

#ifndef __SIGNALS_H__
#define __SIGNALS_H__

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
  size_t field_index;  //!< Index in the sampled data for this filed
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
  double mean;  //!< Mean value
  double sigma;  //!< Standard deviation (constraint)
  bool fixed;  //! Fix the value of the parameter to the mean
};



/**
 * \class Signal
 *
 * A container for signal metadata and PDFs
 */
class Signal {
  public:
    Signal(std::string _name, std::string _title, float _nexpected, float _sigma, std::string _category,
        std::vector<std::string>& hdf5_fields,
        std::vector<size_t>& sample_fields,
        std::vector<Observable>& observables,
        std::vector<Observable>& cuts,
        std::vector<Systematic>& systematics,
        std::vector<std::string>& filenames);
    Signal(std::string _name, std::string _title, float _nexpected, float _sigma, std::string _category,
        std::vector<Observable>& observables,
        std::vector<Systematic>& systematics,
        std::vector<float>& _samples,
        std::vector<int>& _weights);

    std::string name;  //!< string identifier
    std::string title;  //!< histogram title in ROOT-LaTeX format
    std::string category; //!< category like external, cosmogenic, etc for plotting purposes
    double nexpected;  //!< events expected in this fit
    double sigma;  //!< fractional uncertainty
    double efficiency; //!< Fraction of generated events that make it past cuts (not counting the efficiency correction)
    size_t nevents;  //!< number of events in PDF
    pdfz::Eval* histogram;  //!< PDF
};


#endif  // __SIGNALS_H__

