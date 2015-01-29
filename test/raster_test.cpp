#define BOOST_TEST_MODULE lisem raster
#include <boost/test/unit_test.hpp>
#include "raster.h"


BOOST_AUTO_TEST_SUITE(raster)

BOOST_AUTO_TEST_CASE(members)
{
    Raster<int> raster(3, 2);
    raster.cell(0) = -2;
    raster.cell(1) = -1;
    raster.cell(2) = 0;
    raster.cell(3) = 9;
    raster.cell(4) = 1;
    raster.cell(5) = 2;

    BOOST_CHECK_EQUAL(raster.nr_cells(), 6);
    BOOST_CHECK_EQUAL(raster.nr_rows(), 3);
    BOOST_CHECK_EQUAL(raster.nr_cols(), 2);

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
    BOOST_CHECK_EQUAL(raster.nr_cells(), 0u);
    BOOST_CHECK_EQUAL(raster.nr_rows(), 0u);
    BOOST_CHECK_EQUAL(raster.nr_cols(), 0u);
}

BOOST_AUTO_TEST_SUITE_END()
