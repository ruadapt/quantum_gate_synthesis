#define BOOST_TEST_MODULE clifford
#include "../src/clifford.h"
#include <boost/test/included/unit_test.hpp>

namespace c = clifford;

BOOST_AUTO_TEST_CASE(test_creation)
{
    Clifford c1(1, 2, 3, 4);
}