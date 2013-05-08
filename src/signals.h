#ifndef __SIGNALS_H__
#define __SIGNALS_H__

// container for signal metadata and pdfs
struct Signal {
  std::string name;  // string identifier
  std::string title;  // histogram title in ROOT-LaTeX format
  float nexpected;  // events expected in this fit
  float constraint;  // fractional uncertainty
  TH1* histogram;  // PDF (TH1F or TH2F)
};

#endif  // __SIGNALS_H__

