#include <gtest/gtest.h>
#include "hdf5_io.h"

TEST(HDF5IO, ReadWrite1D)
{
	const unsigned int nelements = 5;
	std::vector<float> a(nelements);
	for (unsigned int i=0; i < a.size(); i++)
		a[i] = 2 * i;

	// 1D array
	std::vector<unsigned int> rank(1);
	rank[0] = a.size();

	std::string filename("/tmp/sxmc_test_1d.hdf5");
	std::string dataset("/a");

	ASSERT_TRUE(write_float_vector_hdf5(filename, dataset, a, rank) >= 0);

	///////

	std::vector<float> test_a;
	std::vector<unsigned int> test_rank;

	ASSERT_TRUE(read_float_vector_hdf5(filename, dataset, test_a, test_rank) >= 0);

	ASSERT_EQ((unsigned) 1, test_rank.size());
	EXPECT_EQ(nelements, test_rank[0]);

	ASSERT_EQ(nelements, test_a.size());
	for (unsigned int i=0; i < test_a.size(); i++)
		EXPECT_EQ(2*i, test_a[i]);
}

TEST(HDF5IO, ReadWrite2D)
{
	const unsigned int rows = 3;
	const unsigned int cols = 4;
	const unsigned int nelements = rows * cols;

	std::vector<float> a(nelements);

	for (unsigned int irow=0; irow < rows; irow++)
		for (unsigned int icol=0; icol < cols; icol++)
			a[irow * cols + icol] = 2 * irow + icol;

	// 2D array
	std::vector<unsigned int> rank(2);
	rank[0] = rows;
	rank[1] = cols;

	std::string filename("/tmp/sxmc_test_2d.hdf5");
	std::string dataset("/a");

	ASSERT_TRUE(write_float_vector_hdf5(filename, dataset, a, rank) >= 0);

	///////

	std::vector<float> test_a;
	std::vector<unsigned int> test_rank;

	ASSERT_TRUE(read_float_vector_hdf5(filename, dataset, test_a, test_rank) >= 0);

	ASSERT_EQ((unsigned) 2, test_rank.size());
	EXPECT_EQ(rows, test_rank[0]);
	EXPECT_EQ(cols, test_rank[1]);

	ASSERT_EQ(nelements, test_a.size());

	for (unsigned int irow=0; irow < test_rank[0]; irow++)
		for (unsigned int icol=0; icol < test_rank[1]; icol++)
			EXPECT_EQ(2 * irow + icol, a[irow * test_rank[1] + icol]);
}
