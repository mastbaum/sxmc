#ifndef __OBSERVABLE_H__
#define __OBSERVABLE_H__

/**
 * \file observable.h
 *
 * Fit observables.
*/

#include <string>
#include <vector>

namespace Json {
  class Value;
}

/**
 * \struct Observable
 *
 * A container for observable metadata
*/
struct Observable {
  /** Constructor */
  Observable() {}

  /** Load from JSON configuration */
  Observable(const std::string _name, const Json::Value& config);

  /** Print the observable configuration */
  void print() const;

  std::string name;  //!< Name of the observable
  std::string title;  //!< Title in ROOT LaTeX format, for display
  std::string field;  //!< Name of the field (e.g. "energy")
  std::string units;  //!< Units as string, used in plotting
  std::vector<float> yrange;  //!< y axis range for plots
  size_t field_index;  //!< Index in the sampled data for this field
  size_t bins;  //!< Number of bins
  float lower;  //!< Lower physical bound
  float upper;  //!< Upper physical bound
  bool logscale;  //!< Use log scale y in plots
};

#endif  // __OBSERVABLE_H__

