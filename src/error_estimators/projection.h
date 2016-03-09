/**
 * \file projection.h
 */

#ifndef __ERRORS_PROJECTION_H__
#define __ERRORS_PROJECTION_H__

#include <string>

#include <sxmc/error_estimator.h>
#include <sxmc/projection.h>
#include <sxmc/interval.h>

class LikelihoodSpace;

namespace sxmc {
  namespace errors {

/**
 * 1D projection error estimator.
 *
 * Error estimator which uses a 1D projection of the likelihood space to
 * determine uncertainties. Attempts to put a fraction (1-cl)/2 in both the
 * upper and lower tails. If the mean is less than twice the RMS, push up the
 * upper limit as necessary to assure coverage and report an upper limit.
 */
class Projection : public ErrorEstimator {
  public:
    /**
     * Constructor.
     *
     * \param _lspace The likelihood space samples
     * \param _cl The confidence level
     */
    Projection(LikelihoodSpace* _lspace, float _cl=0.68)
        : ErrorEstimator(_lspace, _cl) {}

    virtual ~Projection() {};

    virtual Interval get_interval(std::string name, float point_estimate);
};

  }  // namespace errors
}  // namespace sxmc

#endif  // __ERRORS_PROJECTION_H__

