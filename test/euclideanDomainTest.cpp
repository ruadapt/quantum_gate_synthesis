#define BOOST_TEST_MODULE euclideanDomain
#include "utils.h"
#include "../euclideanDomain.h"
#include <boost/test/included/unit_test.hpp>

namespace ed = euclidean_domain;

BOOST_AUTO_TEST_CASE(test_rounddiv)
{
    BOOST_CHECK_EQUAL(-2, ed::rounddiv(13, -6));
    BOOST_CHECK_EQUAL(8, ed::rounddiv(100, 12));
    BOOST_CHECK_EQUAL(-14, ed::rounddiv(100, -7));
    BOOST_CHECK_EQUAL(36, ed::rounddiv(-1234, -34));
}