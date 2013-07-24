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
  std::string name;
  std::string title;
  std::string field;
  std::string units;  //!< Units as string, used in plotting
  size_t field_index;
  size_t bins;
  float lower;
  float upper;
};


/**
 * \struct Systematic
 * \brief A container for systematic parameter metadata
 */
struct Systematic {
  std::string name;
  std::string title;
  std::string observable_field;
  std::string truth_field;
  size_t observable_field_index;
  size_t truth_field_index;
  pdfz::Systematic::Type type;
  double mean;
  double sigma;
  bool fixed;
};

#endif  // __SIGNALS_H__

