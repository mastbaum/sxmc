#ifndef __INTERVAL_H__
#define __INTERVAL_H__

/**
 * \file interval.h
 *
 * Confidence/credible intervals.
 */

/** A confidence interval. */
struct Interval {
  /** Constructor, with initialization of parameters to -1. */
  Interval() : one_sided(false), point_estimate(-1), lower(-1), upper(-1),
               cl(-1), coverage(-1) {};

  /** Destructor */
  virtual ~Interval() {}

  /** Get the interval value and uncertainty as a pretty string. */
  std::string str() const;

  bool one_sided;  //!< True if this is a one-sided interval (limit)
  float point_estimate;  //!< Estimate of true value
  float lower;  //!< Lower bound on parameter, for two-sided intervals
  float upper;  //!< Upper bound on parameter
  float cl;  //!< Nominal confidence level
  float coverage;  //!< Actual coverage of this interval
};

#endif  // __INTERVAL_H__

