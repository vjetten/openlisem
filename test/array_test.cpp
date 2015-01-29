#define BOOST_TEST_MODULE lisem array
#include <boost/test/unit_test.hpp>
#include "array.h"


BOOST_AUTO_TEST_SUITE(array)

BOOST_AUTO_TEST_CASE(d1_initializer_list)
{
    Array<int, 1> array{5, 6, 7};

    BOOST_CHECK_EQUAL(array.size(), 3);
    BOOST_CHECK_EQUAL(array.shape().size(), 1);
    BOOST_CHECK_EQUAL(array.shape()[0], 3);
    BOOST_CHECK_EQUAL(array[0], 5);
    BOOST_CHECK_EQUAL(array[1], 6);
    BOOST_CHECK_EQUAL(array[2], 7);
}


BOOST_AUTO_TEST_CASE(d1)
{
    Array<int, 1> array(3);
    array[0] = 5;
    array[1] = 6;
    array[2] = 7;

    BOOST_CHECK_EQUAL(array.size(), 3);
    BOOST_CHECK_EQUAL(array.shape().size(), 1);
    BOOST_CHECK_EQUAL(array.shape()[0], 3);
    BOOST_CHECK_EQUAL(array[0], 5);
    BOOST_CHECK_EQUAL(array[1], 6);
    BOOST_CHECK_EQUAL(array[2], 7);
}


BOOST_AUTO_TEST_CASE(d2)
{
    Array<int, 2> array(3, 2);
    array[0] = -2;
    array[1] = -1;
    array[2] = 0;
    array[3] = 9;
    array[4] = 1;
    array[5] = 2;

    BOOST_CHECK_EQUAL(array.size(), 6);
    BOOST_CHECK_EQUAL(array.shape().size(), 2);
    BOOST_CHECK_EQUAL(array.shape()[0], 3);
    BOOST_CHECK_EQUAL(array.shape()[1], 2);
    BOOST_CHECK_EQUAL(array[0], -2);
    BOOST_CHECK_EQUAL(array[1], -1);
    BOOST_CHECK_EQUAL(array[2], 0);
    BOOST_CHECK_EQUAL(array[3], 9);
    BOOST_CHECK_EQUAL(array[4], 1);
    BOOST_CHECK_EQUAL(array[5], 2);
}


BOOST_AUTO_TEST_CASE(copy)
{
    Array<int, 1> const array{5, 6, 7};
    Array<int, 1> array_copy(array);

    BOOST_CHECK_EQUAL(array.size(), 3);
    BOOST_CHECK_EQUAL(array[0], 5);
    BOOST_CHECK_EQUAL(array[1], 6);
    BOOST_CHECK_EQUAL(array[2], 7);

    BOOST_CHECK_EQUAL(array_copy.size(), 3);
    BOOST_CHECK_EQUAL(array_copy[0], 5);
    BOOST_CHECK_EQUAL(array_copy[1], 6);
    BOOST_CHECK_EQUAL(array_copy[2], 7);
}


BOOST_AUTO_TEST_CASE(move)
{
    Array<int, 1> array{5, 6, 7};
    Array<int, 1> array_move(std::move(array));

    BOOST_CHECK_EQUAL(array.size(), 0);

    BOOST_CHECK_EQUAL(array_move.size(), 3);
    BOOST_CHECK_EQUAL(array_move[0], 5);
    BOOST_CHECK_EQUAL(array_move[1], 6);
    BOOST_CHECK_EQUAL(array_move[2], 7);
}


BOOST_AUTO_TEST_CASE(assign)
{
    Array<int, 1> const array{5, 6, 7};
    Array<int, 1> array_copy = array;

    BOOST_CHECK_EQUAL(array.size(), 3);
    BOOST_CHECK_EQUAL(array[0], 5);
    BOOST_CHECK_EQUAL(array[1], 6);
    BOOST_CHECK_EQUAL(array[2], 7);

    BOOST_CHECK_EQUAL(array_copy.size(), 3);
    BOOST_CHECK_EQUAL(array_copy[0], 5);
    BOOST_CHECK_EQUAL(array_copy[1], 6);
    BOOST_CHECK_EQUAL(array_copy[2], 7);
}


BOOST_AUTO_TEST_CASE(move_assign)
{
    Array<int, 1> array{5, 6, 7};
    Array<int, 1> array_move = std::move(array);

    BOOST_CHECK_EQUAL(array.size(), 0);

    BOOST_CHECK_EQUAL(array_move.size(), 3);
    BOOST_CHECK_EQUAL(array_move[0], 5);
    BOOST_CHECK_EQUAL(array_move[1], 6);
    BOOST_CHECK_EQUAL(array_move[2], 7);
}


BOOST_AUTO_TEST_CASE(empty)
{
    {
        Array<int, 1> array{};
        BOOST_CHECK_EQUAL(array.size(), 0u);
        BOOST_CHECK_EQUAL(array.shape().size(), 1u);
        BOOST_CHECK_EQUAL(array.shape()[0], 0u);
    }

    {
        Array<int, 1> array;
        BOOST_CHECK_EQUAL(array.size(), 0u);
        BOOST_CHECK_EQUAL(array.shape().size(), 1u);
        BOOST_CHECK_EQUAL(array.shape()[0], 0u);
    }

    {
        Array<int, 1> array(0);
        BOOST_CHECK_EQUAL(array.size(), 0u);
        BOOST_CHECK_EQUAL(array.shape().size(), 1u);
        BOOST_CHECK_EQUAL(array.shape()[0], 0u);
    }

    {
        Array<int, 2> array;
        BOOST_CHECK_EQUAL(array.size(), 0u);
        BOOST_CHECK_EQUAL(array.shape().size(), 2u);
        BOOST_CHECK_EQUAL(array.shape()[0], 0u);
        BOOST_CHECK_EQUAL(array.shape()[1], 0u);
    }

    {
        Array<int, 2> array(0, 2);
        BOOST_CHECK_EQUAL(array.size(), 0u);
        BOOST_CHECK_EQUAL(array.shape().size(), 2u);
        BOOST_CHECK_EQUAL(array.shape()[0], 0u);
        BOOST_CHECK_EQUAL(array.shape()[1], 2u);
    }

    {
        Array<int, 2> array(3, 0);
        BOOST_CHECK_EQUAL(array.size(), 0u);
        BOOST_CHECK_EQUAL(array.shape()[0], 3u);
        BOOST_CHECK_EQUAL(array.shape()[1], 0u);
    }
}

BOOST_AUTO_TEST_SUITE_END()
