#ifndef __SYSTEMATIC_H__
#define __SYSTEMATIC_H__

#include <string>
#include <vector>
#include <sxmc/pdfz.h>

namespace Json {
  class Value;
}

/**
 * \struct Systematic
 *
 * A container for systematic parameter metadata
*/
class Systematic {
public:
  Systematic() : means(NULL), sigmas(NULL) {}

  // Load from JSON configuration
  Systematic(const std::string _name, const Json::Value& config);

  virtual ~Systematic() {}

  void print() const;

  std::string name;  //!< Name of the systematic parameter
  std::string title;  //!< Title in ROOT LaTeX format, for display
  std::string observable_field;  //!< Name of the field of the observable
  std::string truth_field;  //! Name of the field of the truth value
  double* means;  //!< Mean values (power series)
  double* sigmas;  //!< Standard deviations
  size_t observable_field_index;  //!< Index of the observable field in the data
  size_t truth_field_index;  //!< Index of the truth field in the data
  size_t npars;  //!< Number of parameters in power series
  std::vector<short> pidx;  //!< Global index for pdfz parameter array offsetting
  pdfz::Systematic::Type type;  //!< The type of systematic
  std::string osc_lut;  //!< Name of oscillation look up table
  bool fixed;  //! Fix the value of the parameter to the mean
};

#endif  // __SYSTEMATIC_H__

