#define BOOST_TEST_MODULE cliffordT
#include "../cliffordT.h"
#include <boost/test/included/unit_test.hpp>

namespace ct = clifford_t;

BOOST_AUTO_TEST_CASE(test_gates_to_string)
{
    List<Gate> gates {X, Y, Z, H, S, T, E, W};
    BOOST_CHECK_EQUAL("XYZHSTEW", ct::to_string(gates));
}