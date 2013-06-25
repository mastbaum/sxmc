#ifndef __CUDA_COMPAT_H__
#define __CUDA_COMPAT_H__

#include <hemi/hemi.h>

#ifndef HEMI_DEV_CODE
// Replacement for atomicAdd on CPU. NOT THREAD-SAFE.
inline unsigned int atomicAdd(unsigned int* address, int val) {
  unsigned int old = *address;
  *address = old + val;
  return old;
}
#endif

#endif  // __CUDA_COMPAT_H__

