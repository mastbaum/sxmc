/**
 * \file hdf5_io.h
 *
 * HDF5 file input and output.
 */

#ifndef __HDF5_IO_H__
#define __HDF5_IO_H__

#include <vector>
#include <string>

namespace sxmc {
  namespace io {

/**
 * Opens an HDF5 file and write the given data to a float
 * dataset.  The rank array is the number of elements in data
 * along each dimension of the array.  This is also recorded
 * in the output file.  The output file is replaced if it already exists.
 *
 * \param filename Name of file to write to
 * \param dataset Name of the HDF5 dataset to write
 * \param data The data to write
 * \param rank The rank of the output data
 * \return Status code, negative in case of failure
 */
int write_float_vector_hdf5(const std::string& filename,
                            const std::string& dataset,
                            const std::vector<float>& data, 
                            const std::vector<unsigned int>& rank);

/**
 * Opens an HDF5 file and reads the given dataset to a float
 * vector.  The rank array is also read from the file.
 * Both data and rank will be resized.
 *
 * \param filename Name of file to write to
 * \param dataset Name of the HDF5 dataset to write
 * \param data The data to write
 * \param rank The rank of the output data
 * \return Status code, negative in case of failure
 */
int read_float_vector_hdf5(const std::string& filename,
                           const std::string& dataset,
                           std::vector<float>& data, 
                           std::vector<unsigned int>& rank);

  }  // namespace io
}  // namespace sxmc

#endif // __HDF5_IO_H__

