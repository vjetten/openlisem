#define BOOST_TEST_MODULE lisem io
#include <boost/test/unit_test.hpp>
#include "CsfMap.h"
#include "error.h"
#include "fixture.h"
#include "io.h"


namespace boost {
namespace test_tools {
namespace tt_detail {

// http://stackoverflow.com/questions/17572583/boost-check-fails-to-compile-operator-for-custom-types
template<>
inline std::ostream& operator<<(
    std::ostream& stream,
    print_helper_t<QString> const& print_helper)
{
    stream << print_helper.m_t.toAscii().constData();
    return stream;
}

} // namespace tt_detail
} // namespace test_tools
} // namespace boost


BOOST_FIXTURE_TEST_SUITE(io, Fixture)

BOOST_AUTO_TEST_CASE(read_raster_does_no_exist)
{
    bool exception_thrown = false;
    ErrorString = "";

    try {
        readRaster("does_not_exist.map");
        exception_thrown = false;
    }
    catch(...) {
        exception_thrown = true;
    }

    BOOST_CHECK(exception_thrown);
    BOOST_CHECK_EQUAL(ErrorString, "Map does_not_exist.map cannot be opened.");
}


BOOST_AUTO_TEST_CASE(read_raster_wrong_format)
{
    bool exception_thrown = false;
    ErrorString = "";

    try {
        readRaster("data/wrong_format.map");
        exception_thrown = false;
    }
    catch(...) {
        exception_thrown = true;
    }

    BOOST_CHECK(exception_thrown);
    BOOST_CHECK_EQUAL(ErrorString,
        "Map data/wrong_format.map cannot be opened.");
}


void test_default_raster(
    std::string const& name)
{
    cTMap raster;
    bool exception_thrown = true;
    ErrorString = "";

    try {
        raster = readRaster(name.c_str());
        exception_thrown = false;
    }
    catch(...) {
        exception_thrown = true;
    }

    BOOST_CHECK_EQUAL(ErrorString, "");
    BOOST_REQUIRE(!exception_thrown);

    BOOST_CHECK_EQUAL(raster.nrRows(), 3);
    BOOST_CHECK_EQUAL(raster.nrCols(), 2);
    BOOST_CHECK_EQUAL(raster.north(), 115.0);
    BOOST_CHECK_EQUAL(raster.west(), -100.0);
    BOOST_CHECK_EQUAL(raster.cellSize(), 5.0);

    BOOST_CHECK(!raster.data.is_mv(0, 0));
    BOOST_CHECK(!raster.data.is_mv(0, 1));
    BOOST_CHECK(!raster.data.is_mv(1, 0));
    BOOST_CHECK( raster.data.is_mv(1, 1));
    BOOST_CHECK(!raster.data.is_mv(2, 0));
    BOOST_CHECK(!raster.data.is_mv(2, 1));

    BOOST_CHECK_EQUAL(raster.data[0][0], -2.0);
    BOOST_CHECK_EQUAL(raster.data[0][1], -1.0);
    BOOST_CHECK_EQUAL(raster.data[1][0],  0.0);
    BOOST_CHECK_EQUAL(raster.data[2][0],  1.0);
    BOOST_CHECK_EQUAL(raster.data[2][1],  2.0);
}


BOOST_AUTO_TEST_CASE(read_raster)
{
    test_default_raster("data/default.map");
    test_default_raster("data/default.img");
}


void test_twice_default_raster(
    std::string const& name)
{
    cTMap raster;
    bool exception_thrown = true;
    ErrorString = "";

    try {
        raster = readRaster(name.c_str());
        exception_thrown = false;
    }
    catch(...) {
        exception_thrown = true;
    }

    BOOST_REQUIRE(!exception_thrown);
    BOOST_CHECK_EQUAL(ErrorString, "");

    BOOST_CHECK_EQUAL(raster.nrRows(), 3);
    BOOST_CHECK_EQUAL(raster.nrCols(), 2);
    BOOST_CHECK_EQUAL(raster.north(), 115.0);
    BOOST_CHECK_EQUAL(raster.west(), -100.0);
    BOOST_CHECK_EQUAL(raster.cellSize(), 5.0);

    BOOST_CHECK(!raster.data.is_mv(0, 0));
    BOOST_CHECK(!raster.data.is_mv(0, 1));
    BOOST_CHECK(!raster.data.is_mv(1, 0));
    BOOST_CHECK( raster.data.is_mv(1, 1));
    BOOST_CHECK(!raster.data.is_mv(2, 0));
    BOOST_CHECK(!raster.data.is_mv(2, 1));

    BOOST_CHECK_EQUAL(raster.data[0][0], -4.0);
    BOOST_CHECK_EQUAL(raster.data[0][1], -2.0);
    BOOST_CHECK_EQUAL(raster.data[1][0],  0.0);
    BOOST_CHECK_EQUAL(raster.data[2][0],  2.0);
    BOOST_CHECK_EQUAL(raster.data[2][1],  4.0);
}


BOOST_AUTO_TEST_CASE(write_map)
{
    auto raster = readRaster("data/default.map");

    raster.data[0][0] *= 2;
    raster.data[0][1] *= 2;
    raster.data[2][0] *= 2;
    raster.data[2][1] *= 2;

    writeRaster(raster, "twice_default.map", "PCRaster");
    test_twice_default_raster("twice_default.map");

    writeRaster(raster, "twice_default.map", "HFA");
    test_twice_default_raster("twice_default.map");
}


BOOST_AUTO_TEST_CASE(write_map_series)
{
    // TODO
}

BOOST_AUTO_TEST_SUITE_END()
