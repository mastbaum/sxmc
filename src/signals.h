#ifndef __SIGNALS_H__
#define __SIGNALS_H__

/**
 * \struct Signal
 * \brief A container for signal metadata and PDFs
 *
 * (Not to be confused with signal.h, which defines UNIX signals.)
 */
struct Signal {
  std::string name;  //!< string identifier
  std::string title;  //!< histogram title in ROOT-LaTeX format
  float nexpected;  //!< events expected in this fit
  float constraint;  //!< fractional uncertainty
  TH1* histogram;  //!< PDF (TH1F or TH2F)
};

#endif  // __SIGNALS_H__

