#ifndef __TTREE_IO_H__
#define __TTREE_IO_H__

/**
 * \file ttree_io.h
 *
 * ROOT TTree input and output.
 */

#include <vector>
#include <string>

namespace sxmc {
  namespace io {

/**
 * Opens a ROOT file and reads the given dataset to a float
 * vector.
 *
 * \param filename - Name of file to read from
 * \param data - The data to read (returned by reference)
 * \param rank - The rank of the input data (returned by reference)
 * \return Status code, negative in case of failure
 */
int read_float_vector_ttree(const std::string& filename,
                            std::vector<float>& data,
                            std::vector<unsigned int>& rank,
                            std::vector<std::string>& ttree_fields);

  }  // namespace io
}  // namespace sxmc

#endif  // __TTREE_IO_H__

