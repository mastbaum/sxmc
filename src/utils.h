/**
 * \file utils.h
 * \brief Collected utility structures and functions
 */

#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>
#include <string>

class TNtuple;

/**
 * A container for a range of values
 */
template <class T>
struct Range {
  T min;  //!< minimum value
  T max;  //!< maximum value
};


/**
 * Get a value from a TNtuple by event ID and field name.
 *
 * \param nt The source TNtuple
 * \param i The event ID
 * \param field The name of the variable to extract
 * \returns The requested value as a float
 */
float get_ntuple_entry(TNtuple* nt, int i, std::string field);


/**
 * Build a correlation matrix for a TNtuple.
 *
 * Creates a matrix with Pearson product-moment correlation coefficients
 * computed between pairs of variables in a TNtuple. The matrix expressed
 * as a vector of length (entries x entries). Only the upper half is set.
 *
 * \param nt The source TNtuple
 * \returns A correlation matrix as a 1D vector
 */
std::vector<float> get_correlation_matrix(TNtuple* nt);


/**
 * \class LinearInterpolator
 * \brief TODO: Document!
 */
class LinearInterpolator {
  public:
    LinearInterpolator() {}

    LinearInterpolator(std::vector<double> _x, std::vector<double> _y)
        : x(_x), y(_y) {}

    void Setup(std::vector<double> _x, std::vector<double> _y) {
      x = _x;
      y = _y;
    }

    double operator()(double _x);

  protected:
    std::vector<double> x;
    std::vector<double> y;
};

#endif  // __UTILS_H__

