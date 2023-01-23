#define BOOST_TEST_MODULE cliffordT
#include "../cliffordT.h"
#include <boost/test/included/unit_test.hpp>

namespace ct = clifford_t;

BOOST_AUTO_TEST_CASE(test_gates_to_string)
{
    List<Gate> gates {X, Y, Z, H, S, T, E, W};
    BOOST_CHECK_EQUAL("XYZHSTEW", ct::to_string(gates));
}

BOOST_AUTO_TEST_CASE(test_twolevel_to_gates)
{
    TwoLevel t1 = make_TL_X(1, 0);
    TwoLevel t2 = make_TL_H(0, 1);
    TwoLevel t3 = make_TL_T(5, 0, 1);
    TwoLevel t4 = make_TL_omega(2, 1);
    BOOST_CHECK_EQUAL(List<Gate>{X}, ct::to_gates(t1));
    BOOST_CHECK_EQUAL(List<Gate>{H}, ct::to_gates(t2));
    BOOST_CHECK_EQUAL((List<Gate>{T, Z}), ct::to_gates(t3));
    BOOST_CHECK_EQUAL(List<Gate>{S}, ct::to_gates(t4));
}

BOOST_AUTO_TEST_CASE(test_twolevels_to_gates)
{
    TwoLevel t1 = make_TL_X(1, 0);
    TwoLevel t2 = make_TL_H(0, 1);
    TwoLevel t3 = make_TL_T(5, 0, 1);
    TwoLevel t4 = make_TL_omega(2, 1);
    BOOST_CHECK_EQUAL((List<Gate>{X, H, T, Z, S}), ct::to_gates(List<TwoLevel>{t1, t2, t3, t4}));
}