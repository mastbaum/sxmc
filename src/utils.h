/**
 * \file utils.h
 * \brief Collected utility structures and functions
 */

#ifndef __UTILS_H__
#define __UTILS_H__

/**
 * \struct Range
 * \brief A container for a range of values
 */
template <class T>
struct Range {
  T min;  //!< minimum value
  T max;  //!< maximum value
};

#endif  // __UTILS_H__

