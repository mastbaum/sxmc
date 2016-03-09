#ifndef __ERROR_ESTIMATOR_H__
#define __ERROR_ESTIMATOR_H__

/**
 * \file error_estimator.h
 *
 * Error (uncertainty) calculation components.
*/

#include <string>

#include <sxmc/interval.h>

class LikelihoodSpace;

/** Types of error estimators. */
typedef enum { ERROR_PROJECTION, ERROR_CONTOUR } ErrorType;


/**
 * Base class for error estimators.
 *
 * Error estimators operate on a LikelihoodSpace to extract parameter
 * uncertainties.
*/
class ErrorEstimator {
  public:
    /**
     * Constructor.
     *
     * \param lspace The likelihood space
     * \param cl Confidence level
    */
    ErrorEstimator(LikelihoodSpace* _lspace, float _cl=0.682)
        : lspace(_lspace), cl(_cl) {}
 
    /** Destructor. */
    virtual ~ErrorEstimator() {}
 
    /**
     * Get the parameter uncertainty.
     *
     * \param name Name of the parameter to get uncertainty for
     * \param point_estimate Estimate of true value
     * \returns An Interval corresponding to the error estimate
    */
    virtual Interval get_interval(std::string name, float point_estimate) = 0;
 
    /** Pretty-print the limit. */
    virtual void print() {};

  protected:
    LikelihoodSpace* lspace;  //!< Likelihood space
    float cl;  //!< Confidence level
};

#endif  // __ERROR_ESTIMATOR_H__

