#ifndef __TESTS_CHI2_H__
#define __TESTS_CHI2_H__

/**
 * \file rms.h
*/

#include <map>
#include <set>
#include <vector>
#include <sxmc/test.h>

class LikelihoodSpace;
class Signal;
struct Source;
struct Observable;
class Systematic;

namespace sxmc {
  namespace tests {

/**
 * \class sxmc::tests::Chi2
 *
 * Compute Chi2/ndf for fits vs. data.
 *
 */
class Chi2 : public Test {
  public:
    /**
     * Constructor.
    */
    Chi2(std::vector<Source>& _sources,
         std::vector<Signal>& _signals,
         std::vector<Systematic>& _systematics,
         std::vector<Observable>& _observables,
         std::set<unsigned>& _datasets,
         std::vector<float>& _samples,
         float _p_threshold=1.5)
      : sources(&_sources), signals(&_signals), systematics(&_systematics),
        observables(&_observables), datasets(&_datasets), samples(&_samples),
        p_threshold(_p_threshold) {}

    /** Destructor */
    virtual ~Chi2() {}

    /** Perform the test. */
    virtual TestResult operator()(const LikelihoodSpace* ls);

  protected:
    std::vector<Source>* sources;
    std::vector<Signal>* signals;
    std::vector<Systematic>* systematics;
    std::vector<Observable>* observables;
    std::set<unsigned>* datasets;
    std::vector<float>* samples;
    double p_threshold;  //!< Maximum Chi2/ndf
};

  }  // namespace tests
}  // namespace sxmc

#endif  // __TESTS_CHI2_H__

