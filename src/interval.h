/**
 * \file interval.h
 *
 * Confidence/credible intervals.
 */

#ifndef __INTERVAL_H__
#define __INTERVAL_H__

/** A confidence interval. */
struct Interval {
  Interval() : one_sided(false), point_estimate(-1), lower(-1), upper(-1),
               cl(-1), coverage(-1) {};

  virtual ~Interval() {}

  /** Get the interval value and uncertainty as a pretty string. */
  std::string str();

  bool one_sided;  //!< True if this is a one-sided interval (limit)
  float point_estimate;  //!< Estimate of true value
  float lower;  //!< Lower bound on parameter, for two-sided intervals
  float upper;  //!< Upper bound on parameter
  float cl;  //!< Nominal confidence level
  float coverage;  //!< Actual coverage of this interval
};

#endif  // __INTERVAL_H__

