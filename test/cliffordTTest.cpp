#define BOOST_TEST_MODULE cliffordT
#include "../src/cliffordT.h"
#include <boost/test/included/unit_test.hpp>

namespace ct = clifford_t;
namespace mat = matrix;

BOOST_AUTO_TEST_CASE(test_gates_to_string)
{
    List<Gate> gates{X, Y, Z, H, S, T, E, W};
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

BOOST_AUTO_TEST_CASE(test_to_gates_syllables_basic)
{
    BOOST_CHECK_EQUAL(List<Gate>{}, ct::to_gates(cons_S_I()));
    BOOST_CHECK_EQUAL(List<Gate>{T}, ct::to_gates(cons_S_T()));
    BOOST_CHECK_EQUAL((List<Gate>{H, T}), ct::to_gates(cons_SApp_HT(cons_S_I())));
    BOOST_CHECK_EQUAL((List<Gate>{T, S, H, T}), ct::to_gates(cons_SApp_SHT(cons_S_T())));
}

BOOST_AUTO_TEST_CASE(test_to_gates_syllables_longer)
{
    List<TrailingSyllable> tail{HT, SHT, SHT, HT, HT, HT, HT, SHT, SHT, HT, HT, HT, HT, HT};
    Syllables ts = Syllables(S_I, tail);
    List<Gate> expected{H, T, H, T, H, T, H, T, H, T, S, H, T, S, H, T, H, T, H, T, H, T, H, T, S, H, T, S, H, T, H, T};
    BOOST_CHECK_EQUAL(expected, ct::to_gates(ts));
}

BOOST_AUTO_TEST_CASE(test_to_gates_clifford)
{
    BOOST_CHECK_EQUAL((List<Gate>{X, S, S, W}), ct::to_gates(Clifford(0, 1, 2, 1)));
}

BOOST_AUTO_TEST_CASE(test_to_gates_normalform)
{
    List<TrailingSyllable> tail{HT, SHT, SHT, HT, HT, HT, HT, SHT, SHT, HT, HT, HT, HT, HT};
    Syllables ts = Syllables(S_I, tail);
    Clifford c = Clifford(0, 1, 2, 1);
    NormalForm n = NormalForm(ts, c);
    List<Gate> expected{H, T, H, T, H, T, H, T, H, T, S, H, T, S, H, T, H, T, H, T, H, T, H, T, S, H, T, S, H, T, H, T, X, S, S, W};
    BOOST_CHECK_EQUAL(expected, ct::to_gates(n));
}

BOOST_AUTO_TEST_CASE(test_normalize)
{
    List<TwoLevel> tls{make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_X(0, 1), make_TL_omega(2, 0), make_TL_omega(5, 1)};
    NormalForm n = ct::normalize(tls);
    List<TrailingSyllable> expected_tail{HT, SHT, SHT, HT, HT, HT, HT, SHT, SHT, HT, HT, HT, HT, HT};
    BOOST_CHECK_EQUAL(S_I, n.ts().lead());
    BOOST_CHECK_EQUAL(expected_tail, n.ts().tail());
    BOOST_CHECK_EQUAL(Clifford(0, 1, 2, 1), n.c());
}

BOOST_AUTO_TEST_CASE(test_synthesis_u2_1)
{
    U2<DOmega> m1 = mat::matrix2x2<DOmega>(1, 0, 0, 1);
    BOOST_CHECK_EQUAL(List<Gate>{}, ct::synthesis_u2(m1));
}

BOOST_AUTO_TEST_CASE(test_synthesis_u2_2)
{
    DOmega a = DOmega(ZDyadic(-1, 1), ZDyadic(-5, 3), ZDyadic(9, 4), ZDyadic(1, 4));
    DOmega b = DOmega(ZDyadic(1, 4), ZDyadic(-1, 3), ZDyadic(1, 3), ZDyadic(-1, 4));
    DOmega c = DOmega(ZDyadic(1, 3), ZDyadic(-1, 3), ZDyadic(1, 4), ZDyadic(1, 4));
    DOmega d = DOmega(ZDyadic(-9, 4), ZDyadic(5, 3), ZDyadic(1, 1), ZDyadic(1, 4));
    U2<DOmega> m = mat::matrix2x2<DOmega>(a, b, c, d);
    List<Gate> expected{H, T, H, T, H, T, H, T, H, T, S, H, T, S, H, T, H, T, H, T, H, T, H, T, S, H, T, S, H, T, H, T, X, S, S, W};
    BOOST_CHECK_EQUAL(expected, ct::synthesis_u2(m));
}