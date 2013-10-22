/**
 * Error estimation things.
 */

#ifndef __ERRORS_H__
#define __ERRORS_H__

#include <string>

class LikelihoodSpace;

/** Types of error estimators. */
typedef enum { ERROR_PROJECTION, ERROR_CONTOUR } ErrorType;


/** A confidence interval. */
struct Interval {
  bool one_sided;  //!< True if this is a one-sided interval (limit)
  float point_estimate;  //!< Estimate of true value
  float lower;  //!< Lower bound on parameter, for two-sided intervals
  float upper;  //!< Upper bound on parameter
  float cl;  //!< Nominal confidence level
  float coverage;  //!< Actual coverage of this interval
};


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
    ErrorEstimator(LikelihoodSpace* _lspace, float _cl=0.68)
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


/**
 * 1D projection error estimator.
 *
 * Error estimator which uses a 1D projection of the likelihood space to
 * determine uncertainties. Attempts to put a fraction (1-cl) in both the
 * upper and lower tails. If including the first bin overcovers, push up the
 * upper limit as necessary to assure coverage and report an upper limit.
 */
class ProjectionError : public ErrorEstimator {
  public:
    ProjectionError(LikelihoodSpace* _lspace, float _cl=0.68)
        : ErrorEstimator(_lspace, _cl) {}

    virtual ~ProjectionError() {};

    virtual Interval get_interval(std::string name, float point_estimate);

  private:
    bool bounded;  //!< Assert a parameter boundary at zero
};


/**
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
class ContourError : public ErrorEstimator {
  public:
    ContourError(LikelihoodSpace* _lspace, float _cl=0.68)
        : ErrorEstimator(_lspace, _cl) {}

    virtual ~ContourError() {};

    virtual Interval get_interval(std::string name, float point_estimate) {
      return Interval();
    }
};

#endif  // __ERRORS_H__

