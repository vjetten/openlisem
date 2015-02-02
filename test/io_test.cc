#define BOOST_TEST_MODULE lisem io
#include <boost/test/unit_test.hpp>
#include "CsfMap.h"
#include "error.h"
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


BOOST_AUTO_TEST_SUITE(io)

BOOST_AUTO_TEST_CASE(load_from_file_does_no_exist)
{
    cTMap raster;
    raster.setPathName("does_not_exist.map");

    bool exception_thrown = false;

    try {
        LoadFromFile(raster);
        exception_thrown = false;
    }
    catch(...) {
        exception_thrown = true;
    }

    BOOST_CHECK(exception_thrown);
    BOOST_CHECK_EQUAL(ErrorString, "Map does_not_exist.map does not exist.");
}


BOOST_AUTO_TEST_CASE(load_from_file_wrong_format)
{
    cTMap raster;
    raster.setPathName("data/wrong_format.map");

    bool exception_thrown = false;

    try {
        LoadFromFile(raster);
        exception_thrown = false;
    }
    catch(...) {
        exception_thrown = true;
    }

    BOOST_CHECK(exception_thrown);
    BOOST_CHECK_EQUAL(ErrorString,
        "Map data/wrong_format.map cannot be opened.");
}


BOOST_AUTO_TEST_CASE(load_from_file)
{
    // Pattern:
    // - Create a cTMap instance.
    // - Set path name.
    // - Read the raster.

    cTMap raster;
    raster.setPathName("data/default.map");

    bool success = false;
    bool exception_thrown = true;

    try {
        success = LoadFromFile(raster);
        exception_thrown = false;
    }
    catch(...) {
        BOOST_CHECK_EQUAL(ErrorString, "");
        exception_thrown = true;
    }

    BOOST_REQUIRE(!exception_thrown);
    BOOST_REQUIRE(success);

    BOOST_CHECK_EQUAL(raster.nrRows(), 3);
    BOOST_CHECK_EQUAL(raster.nrCols(), 2);
    BOOST_CHECK_EQUAL(raster.north(), 115.0);
    BOOST_CHECK_EQUAL(raster.west(), -100.0);
    BOOST_CHECK_EQUAL(raster.cellSize(), 5.0);

    BOOST_CHECK(!raster.Data.is_mv(0, 0));
    BOOST_CHECK(!raster.Data.is_mv(0, 1));
    BOOST_CHECK(!raster.Data.is_mv(1, 0));
    BOOST_CHECK( raster.Data.is_mv(1, 1));
    BOOST_CHECK(!raster.Data.is_mv(2, 0));
    BOOST_CHECK(!raster.Data.is_mv(2, 1));

    BOOST_CHECK_EQUAL(raster.Data[0][0], -2.0);
    BOOST_CHECK_EQUAL(raster.Data[0][1], -1.0);
    BOOST_CHECK_EQUAL(raster.Data[1][0],  0.0);
    BOOST_CHECK_EQUAL(raster.Data[2][0],  1.0);
    BOOST_CHECK_EQUAL(raster.Data[2][1],  2.0);
}


BOOST_AUTO_TEST_CASE(write_map)
{
    // Pattern:
    // - Create a cTMap instance.
    // - Set path name.
    // - Write the raster.
    cTMap raster;
    raster.setPathName("data/default.map");
    LoadFromFile(raster);

    raster.Data[0][0] *= 2;
    raster.Data[0][1] *= 2;
    raster.Data[2][0] *= 2;
    raster.Data[2][1] *= 2;

    WriteMap(raster, "twice_default.map");

    raster.setPathName("twice_default.map");

    bool success = false;
    bool exception_thrown = true;

    try {
        success = LoadFromFile(raster);
        exception_thrown = false;
    }
    catch(...) {
        BOOST_CHECK_EQUAL(ErrorString, "");
        exception_thrown = true;
    }

    BOOST_REQUIRE(!exception_thrown);
    BOOST_REQUIRE(success);

    BOOST_CHECK_EQUAL(raster.nrRows(), 3);
    BOOST_CHECK_EQUAL(raster.nrCols(), 2);
    BOOST_CHECK_EQUAL(raster.north(), 115.0);
    BOOST_CHECK_EQUAL(raster.west(), -100.0);
    BOOST_CHECK_EQUAL(raster.cellSize(), 5.0);

    BOOST_CHECK(!raster.Data.is_mv(0, 0));
    BOOST_CHECK(!raster.Data.is_mv(0, 1));
    BOOST_CHECK(!raster.Data.is_mv(1, 0));
    BOOST_CHECK( raster.Data.is_mv(1, 1));
    BOOST_CHECK(!raster.Data.is_mv(2, 0));
    BOOST_CHECK(!raster.Data.is_mv(2, 1));

    BOOST_CHECK_EQUAL(raster.Data[0][0], -4.0);
    BOOST_CHECK_EQUAL(raster.Data[0][1], -2.0);
    BOOST_CHECK_EQUAL(raster.Data[1][0],  0.0);
    BOOST_CHECK_EQUAL(raster.Data[2][0],  2.0);
    BOOST_CHECK_EQUAL(raster.Data[2][1],  4.0);
}


BOOST_AUTO_TEST_CASE(write_map_series)
{
    // TODO
}

BOOST_AUTO_TEST_SUITE_END()
