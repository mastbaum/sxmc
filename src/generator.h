/**
 * \class FakeDataGenerator
 * \brief Generate fake data sets
 *
 * Make fake data set Ntuples by sampling PDFs according to given
 * normalizations.
 */

#ifndef __GENERATOR_H__
#define __GENERATOR_H__

#include "utils.h"
#include "signals.h"

class FakeDataGenerator {
  public:
    /**
     * Constructor
     *
     * \param signals vector of signals whose PDFs to samples
     * \param _e_range Energy range in which to generate events
     * \param _r_range Radius range in which to generate events
     */
    FakeDataGenerator(std::vector<Signal> signals, Range<float> _e_range,
                      Range<float> _r_range)
      : e_range(_e_range), r_range(_r_range) {
      for (auto it=signals.cbegin(); it!=signals.cend(); it++) {
        this->pdfs.push_back(it->histogram);
      }
    }

    virtual ~FakeDataGenerator() {}

    /**
     * Create a dataset (an ntuple with fields "r:e")
     *
     * If histograms are TH2s, both fields are filled; if TH1, r is set to 0.
     */
    TNtuple* operator()(float* norms, bool poisson=true);

  protected:
    std::vector<TH1*> pdfs;  //!< The set of PDF histograms
    Range<float> e_range;  //!< Energy range in which to generate events
    Range<float> r_range;  //!< Radius range in which to generate events
};

#endif  // __GENERATOR_H__

