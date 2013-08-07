#ifndef __SIGNALS_H__
#define __SIGNALS_H__

#include "pdfz.h"

/**
 * \struct Signal
 * \brief A container for signal metadata and PDFs
 *
 * (Not to be confused with signal.h, which defines UNIX signals.)
 */
struct Signal {
  std::string name;  //!< string identifier
  std::string title;  //!< histogram title in ROOT-LaTeX format
  double nexpected;  //!< events expected in this fit
  double sigma;  //!< fractional uncertainty
  size_t nevents;  //!< number of events in PDF
  pdfz::Eval* histogram;  //!< PDF
};


/**
 * \struct Observable
 * \brief A container for observable metadata
 */
struct Observable {
  std::string name;  //!< Name of the observable
  std::string title;  //!< Title in ROOT LaTeX format, for display
  std::string field;  //!< Name of the field (e.g. "energy")
  std::string units;  //!< Units as string, used in plotting
  size_t field_index;  //!< Index of data field in the HDF5 table
  size_t bins;  //!< Number of bins
  float lower;  //!< Lower physical bound
  float upper;  //!< Upper physical bound
  float exclude_min;  //!< Minimum of excluded window
  float exclude_max;  //!< Maximum of excluded window
  bool exclude;  //!< Exclude a window inside the fit range
};


/**
 * \struct Systematic
 * \brief A container for systematic parameter metadata
 */
struct Systematic {
  std::string name;  //!< Name of the systematic parameter
  std::string title;  //!< Title in ROOT LaTeX format, for display
  std::string observable_field;  //!< Name of the field of the observable
  std::string truth_field;  //! Name of the field of the truth value
  size_t observable_field_index;  //!< Index of the observable field in HDF5
  size_t truth_field_index;  //!< Index of the truth field in HDF5
  pdfz::Systematic::Type type;  //!< The type of systematic
  double mean;  //!< Mean value
  double sigma;  //!< Standard deviation (constraint)
  bool fixed;  //! Fix the value of the parameter to the mean
};

#endif  // __SIGNALS_H__

