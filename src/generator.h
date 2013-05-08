#ifndef __GENERATOR_H__
#define __GENERATOR_H__

#include "utils.h"
#include "signals.h"

// create fake datasets by sampling histograms
class FakeDataGenerator {
  public:
    FakeDataGenerator(std::vector<Signal> signals, Range<float> _e_range,
                      Range<float> _r_range)
      : e_range(_e_range), r_range(_r_range) {
      for (auto it=signals.cbegin(); it!=signals.cend(); it++) {
        this->pdfs.push_back(it->histogram);
      }
    }
    virtual ~FakeDataGenerator() {}

    // create a dataset: an ntuple with fields "r:e"
    // if histograms are TH2, both are filled; if TH1, r is set to 0.
    TNtuple* operator()(float* norms, bool poisson=true);

  protected:
    std::vector<TH1*> pdfs;
    Range<float> e_range;
    Range<float> r_range;
};

#endif  // __GENERATOR_H__

