/**
 * \file ttree_io.h
 *
 * root TTree input and output.
 */

#ifndef __TTREE_IO_H__
#define __TTREE_IO_H__

#include <vector>
#include <string>

// because c++ is the devil
struct my_bool
{
      bool the_bool;
};
 

/**
 * Opens a root file and reads the given dataset to a float
 * vector.
 *
 * \param filename Name of file to read from
 * \param data The data to read
 * \param rank The rank of the input data
 * \return Status code, negative in case of failure
 */

int read_float_vector_ttree(const std::string &filename,
                            std::vector<float> &data, 
                            std::vector<unsigned int> &rank,
                            std::vector<std::string> &ttree_fields);

#endif // __TTREE_IO_H__

