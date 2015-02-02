#define BOOST_TEST_MODULE lisem raster
#include <boost/test/unit_test.hpp>
#include "raster.h"


BOOST_AUTO_TEST_SUITE(raster)

BOOST_AUTO_TEST_CASE(members)
{
    size_t nr_rows{3};
    size_t nr_cols{2};
    double north{100.0};
    double west{-100.0};
    double cell_size{50.0};

    Raster<int> raster(nr_rows, nr_cols, north, west, cell_size);
    raster.cell(0) = -2;
    raster.cell(1) = -1;
    raster.cell(2) = 0;
    raster.cell(3) = 9;
    raster.cell(4) = 1;
    raster.cell(5) = 2;

    BOOST_CHECK_EQUAL(raster.nr_rows(), nr_rows);
    BOOST_CHECK_EQUAL(raster.nr_cols(), nr_cols);
    BOOST_CHECK_EQUAL(raster.nr_cells(), nr_rows * nr_cols);
    BOOST_CHECK_EQUAL(raster.north(), north);
    BOOST_CHECK_EQUAL(raster.west(), west);
    BOOST_CHECK_EQUAL(raster.cell_size(), cell_size);

    BOOST_CHECK_EQUAL(raster[0][0], -2);
    BOOST_CHECK_EQUAL(raster[0][1], -1);
    BOOST_CHECK_EQUAL(raster[1][0], 0);
    BOOST_CHECK_EQUAL(raster[1][1], 9);
    BOOST_CHECK_EQUAL(raster[2][0], 1);
    BOOST_CHECK_EQUAL(raster[2][1], 2);

    BOOST_CHECK_EQUAL(raster.cell(0, 0), -2);
    BOOST_CHECK_EQUAL(raster.cell(0, 1), -1);
    BOOST_CHECK_EQUAL(raster.cell(1, 0), 0);
    BOOST_CHECK_EQUAL(raster.cell(1, 1), 9);
    BOOST_CHECK_EQUAL(raster.cell(2, 0), 1);
    BOOST_CHECK_EQUAL(raster.cell(2, 1), 2);

    BOOST_CHECK_EQUAL(raster.cell(0), -2);
    BOOST_CHECK_EQUAL(raster.cell(1), -1);
    BOOST_CHECK_EQUAL(raster.cell(2), 0);
    BOOST_CHECK_EQUAL(raster.cell(3), 9);
    BOOST_CHECK_EQUAL(raster.cell(4), 1);
    BOOST_CHECK_EQUAL(raster.cell(5), 2);
}


BOOST_AUTO_TEST_CASE(default_constructor)
{
    Raster<int> raster;
    BOOST_CHECK_EQUAL(raster.nr_rows(), 0u);
    BOOST_CHECK_EQUAL(raster.nr_cols(), 0u);
    BOOST_CHECK_EQUAL(raster.nr_cells(), 0u);
    BOOST_CHECK_EQUAL(raster.north(), 0u);
    BOOST_CHECK_EQUAL(raster.west(), 0u);
    BOOST_CHECK_EQUAL(raster.cell_size(), 0u);
}

BOOST_AUTO_TEST_SUITE_END()
