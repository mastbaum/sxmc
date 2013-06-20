#include "hdf5_io.h"
#include <hdf5.h>
#include <hdf5_hl.h>

herr_t write_float_vector_hdf5(const std::string &filename,
                               const std::string &dataset,
                               const std::vector<float> &data, 
                               const std::vector<unsigned int> &rank)
{
    hid_t       file_id;
    herr_t      status;

    std::vector<hsize_t> dims(rank.size());
    for (unsigned int i=0; i < rank.size(); i++)
        dims[i] = rank[i];

    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    if (file_id < 0)
        return file_id; // Fail

    status = H5LTmake_dataset_float(file_id, dataset.c_str(), dims.size(), &dims.front(), &data.front());

    // Do regardless in case above didn't work
    H5Fclose(file_id);

    return status;
}


herr_t read_float_vector_hdf5(const std::string &filename,
                              const std::string &dataset,
                              std::vector<float> &data, 
                              std::vector<unsigned int> &rank)
{
    hid_t       file_id;
    herr_t      status;
    int ndims;
    std::vector<hsize_t> dims;
    unsigned int nelements = 1;

    // Open file
    file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    if (file_id < 0)
        return file_id; // Fail

    // Get number of dimensions in array
    status = H5LTget_dataset_ndims(file_id, dataset.c_str(), &ndims);
    if (status < 0)
        goto close_file;  // Yeah, I went there.

    // Get rank of each dimension
    dims.resize(ndims);
    status = H5LTget_dataset_info(file_id, dataset.c_str(), &dims.front(), NULL, NULL);
    if (status < 0)
        goto close_file; // No really, Linus said it was OK.

    rank.resize(ndims);
    for (unsigned int i=0; i < rank.size(); i++) {
        rank[i] = dims[i];
        nelements *= dims[i];
    }

    // Read the data
    data.resize(nelements);
    status = H5LTread_dataset_float(file_id, dataset.c_str(), &data.front());

close_file:
    H5Fclose(file_id);
    return status;
}
