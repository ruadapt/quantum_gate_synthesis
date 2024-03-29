#define BOOST_TEST_MODULE multiQubitSynthesis
#include "../src/multiQubitSynthesis.h"
#include <boost/test/included/unit_test.hpp>

namespace mqs = multi_qubit_synthesis;
namespace mat = matrix;

BOOST_AUTO_TEST_CASE(test_invert_twolevel)
{
    TwoLevel t1 = make_TL_X(5, 6);
    TwoLevel t2 = make_TL_H(7, 8);
    TwoLevel t3 = make_TL_T(9, 10, 11);
    TwoLevel t4 = make_TL_omega(12, 13);
    BOOST_CHECK_EQUAL(make_TL_X(5, 6), mqs::invert_twolevel(t1));
    BOOST_CHECK_EQUAL(make_TL_H(7, 8), mqs::invert_twolevel(t2));
    BOOST_CHECK_EQUAL(make_TL_T(-9, 10, 11), mqs::invert_twolevel(t3));
    BOOST_CHECK_EQUAL(make_TL_omega(-12, 13), mqs::invert_twolevel(t4));
}

BOOST_AUTO_TEST_CASE(test_twolevel_matrix)
{
    Matrix<int, 3, 3> m = mqs::twolevel_matrix<int, 3>(Pair<int>{5, 6}, Pair<int>{7, 8}, 1, 2);
    BOOST_CHECK_EQUAL(1, m(0, 0));
    BOOST_CHECK_EQUAL(0, m(0, 1));
    BOOST_CHECK_EQUAL(0, m(0, 2));
    BOOST_CHECK_EQUAL(0, m(1, 0));
    BOOST_CHECK_EQUAL(5, m(1, 1));
    BOOST_CHECK_EQUAL(6, m(1, 2));
    BOOST_CHECK_EQUAL(0, m(2, 0));
    BOOST_CHECK_EQUAL(7, m(2, 1));
    BOOST_CHECK_EQUAL(8, m(2, 2));
}

BOOST_AUTO_TEST_CASE(test_onelevel_matrix)
{
    Matrix<int, 3, 3> m = mqs::onelevel_matrix<int, 3>(5, 2);
    BOOST_CHECK_EQUAL(1, m(0, 0));
    BOOST_CHECK_EQUAL(0, m(0, 1));
    BOOST_CHECK_EQUAL(0, m(0, 2));
    BOOST_CHECK_EQUAL(0, m(1, 0));
    BOOST_CHECK_EQUAL(1, m(1, 1));
    BOOST_CHECK_EQUAL(0, m(1, 2));
    BOOST_CHECK_EQUAL(0, m(2, 0));
    BOOST_CHECK_EQUAL(0, m(2, 1));
    BOOST_CHECK_EQUAL(5, m(2, 2));
}

BOOST_AUTO_TEST_CASE(test_matrix_of_twolevel)
{
    TwoLevel t1 = make_TL_X(1, 2);
    TwoLevel t2 = make_TL_H(2, 2);
    TwoLevel t3 = make_TL_T(1, 2, 2);
    TwoLevel t4 = make_TL_omega(2, 1);

    Matrix<DOmega, 2, 2> m1 = mqs::matrix_of_twolevel<DOmega, 2>(t1);
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 1), m1(0, 0));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m1(0, 1));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m1(1, 0));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m1(1, 1));

    Matrix<DOmega, 2, 2> m2 = mqs::matrix_of_twolevel<DOmega, 2>(t2);
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 1), m2(0, 0));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m2(0, 1));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m2(1, 0));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 1), m2(1, 1));

    Matrix<DOmega, 2, 2> m3 = mqs::matrix_of_twolevel<DOmega, 2>(t3);
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 1), m3(0, 0));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m3(0, 1));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m3(1, 0));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 1), m3(1, 1));

    Matrix<DOmega, 2, 2> m4 = mqs::matrix_of_twolevel<DOmega, 2>(t4);
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 1), m4(0, 0));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m4(0, 1));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m4(1, 0));
    BOOST_CHECK_EQUAL(DOmega(0, 1, 0, 0), m4(1, 1));
}

BOOST_AUTO_TEST_CASE(test_matrix_of_twolevels)
{
    TwoLevel t1 = make_TL_X(1, 5);
    TwoLevel t2 = make_TL_H(2, 2);
    TwoLevel t3 = make_TL_T(1, 2, 2);
    TwoLevel t4 = make_TL_omega(2, 1);

    Matrix<DOmega, 3, 3> m = mqs::matrix_of_twolevels<DOmega, 3>(List<TwoLevel>{t1, t2, t3, t4});
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 1), m(0, 0));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m(0, 1));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m(0, 2));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m(1, 0));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m(1, 1));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m(1, 2));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m(2, 0));
    BOOST_CHECK_EQUAL(DOmega(0, 0, 0, 0), m(2, 1));
    BOOST_CHECK_EQUAL(DOmega(ZDyadic(-1, 1), 0, ZDyadic(1, 1), 0), m(2, 2));
}

BOOST_AUTO_TEST_CASE(test_reduce_column)
{
    Matrix<DOmega, 2, 1> v;
    v(0, 0) = DOmega(ZDyadic(1, 2), 0, ZDyadic(1, 3), 0);
    v(1, 0) = DOmega(ZDyadic(5, 3), 0, ZDyadic(5, 3), ZDyadic(3, 3));
    Index i = 0;
    List<TwoLevel> gs = mqs::reduce_column(v, i);
    List<TwoLevel> expected{make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_omega(0, 0)};
    BOOST_CHECK_EQUAL(expected, gs);
}

BOOST_AUTO_TEST_CASE(test_synthesis_nqubit)
{
    U2<DOmega> m = mat::matrix2x2(1, 0, 0, 1);
    List<TwoLevel> gs = mqs::synthesis_nqubit(m);
    List<TwoLevel> expected{make_TL_omega(0, 0), make_TL_omega(0, 1)};
    BOOST_CHECK_EQUAL(expected, gs);
}

BOOST_AUTO_TEST_CASE(test_synthesis_nqubit2)
{
    DOmega u = DOmega(ZDyadic(1, 2), 0, ZDyadic(1, 3), 0);
    DOmega t = DOmega(ZDyadic(5, 3), 0, ZDyadic(5, 3), ZDyadic(3, 3));
    U2<DOmega> m = mat::matrix2x2(u, -t.adj(), t, u.adj());
    List<TwoLevel> gs = mqs::synthesis_nqubit(m);
    List<TwoLevel> expected{make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_omega(0, 0), make_TL_omega(1, 1)};
    BOOST_CHECK_EQUAL(expected, gs);
}

BOOST_AUTO_TEST_CASE(test_synthesis_nqubit3)
{
    DOmega a = DOmega(ZDyadic(-1, 1), ZDyadic(-5, 3), ZDyadic(9, 4), ZDyadic(1, 4));
    DOmega b = DOmega(ZDyadic(1, 4), ZDyadic(-1, 3), ZDyadic(1, 3), ZDyadic(-1, 4));
    DOmega c = DOmega(ZDyadic(1, 3), ZDyadic(-1, 3), ZDyadic(1, 4), ZDyadic(1, 4));
    DOmega d = DOmega(ZDyadic(-9, 4), ZDyadic(5, 3), ZDyadic(1, 1), ZDyadic(1, 4));
    U2<DOmega> m = mat::matrix2x2<DOmega>(a, b, c, d);
    List<TwoLevel> gs = mqs::synthesis_nqubit(m);
    List<TwoLevel> expected{make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_T(3, 0, 1), make_TL_H(0, 1), make_TL_T(1, 0, 1), make_TL_H(0, 1), make_TL_X(0, 1), make_TL_omega(2, 0), make_TL_omega(5, 1)};
    BOOST_CHECK_EQUAL(expected, gs);
}