#define BOOST_TEST_MODULE matrix
#include "../matrix.h"
#include "../ring.h"
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/multiprecision/number.hpp>

BOOST_AUTO_TEST_CASE(test_indexing)
{
    Vector<QRootTwo, 7> v;
    v(1) = QRootTwo(5_mpq / 6, 7_mpq / 8);
    assert(QRootTwo(5_mpq / 6, 7_mpq / 8) == v(1));
    assert(QRootTwo(0) == v(2));

    Matrix<ZRootTwo, 4, 4> m;
    m(1, 1) = ZRootTwo(4, 6);
    assert(ZRootTwo(4, 6) == m(1, 1));
    assert(ZRootTwo(0) == m(1, 2));
}

BOOST_AUTO_TEST_CASE(test_matrix_addition_int)
{
    Matrix<int, 2, 2> m;
    m(0, 0) = 5;
    m(0, 1) = 7;
    m(1, 0) = 2;
    m(1, 1) = 3;
    Matrix<int, 2, 2> m2;
    m2(0, 0) = 3;
    m2(0, 1) = 3;
    m2(1, 0) = 2;
    m2(1, 1) = 1;
    Matrix<int, 2, 2> sum = m + m2;
    assert(8 == sum(0, 0));
    assert(10 == sum(0, 1));
    assert(4 == sum(1, 0));
    assert(4 == sum(1, 1));
}

BOOST_AUTO_TEST_CASE(test_matrix_addition_Integer)
{
    Matrix<Integer, 2, 2> m;
    m(0, 0) = 5;
    m(0, 1) = 7;
    m(1, 0) = 2;
    m(1, 1) = 3;
    Matrix<Integer, 2, 2> m2;
    m2(0, 0) = 3;
    m2(0, 1) = 3;
    m2(1, 0) = 2;
    m2(1, 1) = 1;
    Matrix<Integer, 2, 2> sum = m + m2;
    assert(8_mpz == sum(0, 0));
    assert(10_mpz == sum(0, 1));
    assert(4_mpz == sum(1, 0));
    assert(4_mpz == sum(1, 1));
}

BOOST_AUTO_TEST_CASE(test_matrix_addition_QRootTwo)
{
    Matrix<QRootTwo, 2, 2> m;
    m(0, 0) = 5;
    m(0, 1) = 7;
    m(1, 0) = 2;
    m(1, 1) = 3;
    Matrix<QRootTwo, 2, 2> m2;
    m2(0, 0) = 3;
    m2(0, 1) = 3;
    m2(1, 0) = 2;
    m2(1, 1) = 1;
    Matrix<QRootTwo, 2, 2> sum = m + m2;
    assert(QRootTwo(8) == sum(0, 0));
    assert(QRootTwo(10) == sum(0, 1));
    assert(QRootTwo(4) == sum(1, 0));
    assert(QRootTwo(4) == sum(1, 1));
}

BOOST_AUTO_TEST_CASE(test_matrix_prod_QRootTwo)
{
    Matrix<QRootTwo, 2, 2> m;
    m(0, 0) = 5;
    m(0, 1) = 7;
    m(1, 0) = 2;
    m(1, 1) = 3;
    Matrix<QRootTwo, 2, 2> m2;
    m2(0, 0) = 3;
    m2(0, 1) = 3;
    m2(1, 0) = 2;
    m2(1, 1) = QRootTwo(1_mpq / 2, 2_mpq / 3);
    Matrix<QRootTwo, 2, 2> p = prod(m, m2);
    assert(QRootTwo(29) == p(0, 0));
    assert(QRootTwo(37_mpq / 2, 14_mpq / 3) == p(0, 1));
    assert(QRootTwo(12) == p(1, 0));
    assert(QRootTwo(15_mpq / 2, 2) == p(1, 1));
}

BOOST_AUTO_TEST_CASE(test_matrix_of_function)
{
    auto f = [](size_t i, size_t j) -> int
    {
        return int(i + j);
    };
    Matrix<int, 3, 4> m = matrix_of_function<int, 3, 4>(f);
    BOOST_CHECK_EQUAL(0, m(0, 0));
    BOOST_CHECK_EQUAL(1, m(0, 1));
    BOOST_CHECK_EQUAL(2, m(0, 2));
    BOOST_CHECK_EQUAL(3, m(0, 3));
    BOOST_CHECK_EQUAL(1, m(1, 0));
    BOOST_CHECK_EQUAL(2, m(1, 1));
    BOOST_CHECK_EQUAL(3, m(1, 2));
    BOOST_CHECK_EQUAL(4, m(1, 3));
    BOOST_CHECK_EQUAL(2, m(2, 0));
    BOOST_CHECK_EQUAL(3, m(2, 1));
    BOOST_CHECK_EQUAL(4, m(2, 2));
    BOOST_CHECK_EQUAL(5, m(2, 3));
}