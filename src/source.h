#ifndef __SOURCE_H__
#define __SOURCE_H__

/**
 * \file source.h
 *
 * Fit sources.
*/

#include <string>

namespace Json {
  class Value;
}

/**
 * \struct Source
 *
 * A container for a signal source
*/
struct Source {
  /** Default constructor */
  Source() {}

  /**
   * Constructor.
   *
   * Specify source parameters explicitly.
   *
   * \param _name - Name of the source
   * \param _index - Index in the list of sources
   * \param _mean - Mean value for the source rate
   * \param _sigma - Gaussian constraint on the rate
   * \param _fixed - Fix the rate in the fit
  */
  Source(const std::string& _name, size_t _index,
         float _mean, float _sigma, bool _fixed)
      : name(_name), index(_index),
        mean(_mean), sigma(_sigma), fixed(_fixed) {}

  /** Constructor.
   *
   * Load parameters from a JSON configuration.
   *
   * \param _name - Name of the source
   * \param params - JSON object containing parameters
  */
  Source(const std::string& _name, const Json::Value& params);

  /** Print the source configuration. */
  void print() const;

  std::string name;  //!< String identifier
  size_t index;  //!< Index in the list of sources
  float mean;  //!< Mean expectation (scaling, 1.0 is nominal)
  float sigma;  //!< Gaussian constraint (fractional)
  bool fixed;  //!< Rate is fixed
};

#endif  // __SOURCE_H__

