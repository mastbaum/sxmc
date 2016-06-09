#ifndef __SOURCE_H__
#define __SOURCE_H__

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
  Source() {}

  Source(const std::string& _name, size_t _index,
         float _mean, float _sigma, bool _fixed)
      : name(_name), index(_index),
        mean(_mean), sigma(_sigma), fixed(_fixed) {}

  // Load from JSON configuration
  Source(const std::string& _name, const Json::Value& params);

  void print() const;

  std::string name;  //!< String identifier
  size_t index;  //!< Index in the list of sources
  float mean;  //!< Mean expectation (scaling, 1.0 is nominal)
  float sigma;  //!< Gaussian constraint (fractional)
  bool fixed;  //!< Rate is fixed
};

#endif  // __SOURCE_H__

