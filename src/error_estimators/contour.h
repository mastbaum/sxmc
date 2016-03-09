#ifndef __ERRORS_CONTOUR_H__
#define __ERRORS_CONTOUR_H__

/**
 * \file contour.h
*/

#include <string>

#include <sxmc/contour.h>
#include <sxmc/error_estimator.h>
#include <sxmc/interval.h>

class TNtuple;
class LikelihoodSpace;

namespace sxmc {
  namespace errors {

/**
 * \class sxmc::errors::Contour
 *
 * Countour/profile likelihood error estimator.
 *
 * Error estimator which uses the boundary of an n-dimensional contour in the
 * full likelihood space to estimate uncertainties, which is equivalent to the
 * profile likelihood construction.
 *
 * In practice, this means finding the likelihood surface which contains points
 * within N units of the maximum and taking the minimum and maximum in each
 * dimension, where N is half of the quantile of the chi squared distribution
 * corresponding to the desired confidence level.
 */
class Contour : public ErrorEstimator {
  public:
    /**
     * Constructor
     *
     * \param _lspace - The LikelihoodSpace on which to operate
     * \param _cl - The confidence level for error estimation
    */
    Contour(LikelihoodSpace* _lspace, float _cl=0.68);

    /** Destructor */
    virtual ~Contour();

    /**
     * Get the interval for a given parameter.
     *
     * \param name - The name of the parameter
     * \param point_estimate - The best-fit point
     * \returns An Interval with the error estimate
    */
    virtual Interval get_interval(std::string name, float point_estimate);

  protected:
    TNtuple* contour_points;  //!< Likelihood space samples within the contour
};

  }  // namespace errors
}  // namespace sxmc

#endif  // __ERRORS_CONTOUR_H__

