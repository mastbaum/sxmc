#ifndef __TESTS_DRIFT_H__
#define __TESTS_DRIFT_H__

/**
 * \file rms.h
*/

#include <sxmc/test.h>

class LikelihoodSpace;

namespace sxmc {
  namespace tests {

/**
 * \class sxmc::test::Drift
 *
 * A test that the mean and RMS stays constant.
 *
 * In a well-mixed chain, the mean should not wander far beyond the RMS.
 * For this test, we slice the samples into N chunks and make sure the
 * parameter distribution is consistent over time by calculting Chi^2
 * p-values for each slice against each other slice.
 */
class Drift : public Test {
  public:
    /**
     * Constructor
     *
     * \param _ks_threshold - Minimum probability
     * \param _slices - Number of time slices
    */
    Drift(unsigned _slices, float _p_threshold=0.005)
      : slices(_slices), p_threshold(_p_threshold) {}

    /** Destructor */
    virtual ~Drift() {}

    /** Perform the test. */
    virtual TestResult operator()(const LikelihoodSpace* ls);

  protected:
    unsigned slices;  //!< Number of time slices
    double p_threshold;  //!< Minimum probability
};

  }  // namespace tests
}  // namespace sxmc

#endif  // __TESTS_DRIFT_H__

