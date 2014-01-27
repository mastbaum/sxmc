/**
 * \file cuda_compat.h
 *
 * Definitions required for compatibility between CPU and GPU targets.
 */
#ifndef __CUDA_COMPAT_H__
#define __CUDA_COMPAT_H__

#include <hemi/hemi.h>

#ifndef HEMI_DEV_CODE
/**
 * Replacement for atomicAdd on CPU. NOT THREAD-SAFE.
 *
 * \param address Address of the integer to add to
 * \param val Value to add to the integer at address
 * \returns The value stored at the address before the addition operation
 */
inline unsigned int atomicAdd(unsigned int* address, int val) {
  unsigned int old = *address;
  *address = old + val;
  return old;
}
#endif

#endif  // __CUDA_COMPAT_H__

