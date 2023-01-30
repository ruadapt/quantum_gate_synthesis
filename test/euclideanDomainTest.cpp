#define BOOST_TEST_MODULE euclideanDomain
#include "comparisons.h"
#include "../src/euclideanDomain.h"
#include <boost/test/included/unit_test.hpp>

namespace ed = euclidean_domain;

BOOST_AUTO_TEST_CASE(test_rounddiv)
{
    BOOST_CHECK_EQUAL(-2, ed::rounddiv(13, -6));
    BOOST_CHECK_EQUAL(8, ed::rounddiv(100, 12));
    BOOST_CHECK_EQUAL(-14, ed::rounddiv(100, -7));
    BOOST_CHECK_EQUAL(36, ed::rounddiv(-1234, -34));
}

BOOST_AUTO_TEST_CASE(test_divMod)
{
    std::tuple<Integer, Integer> t1 = ed::divMod(100, 12);
    BOOST_CHECK_EQUAL(8, std::get<0>(t1));
    BOOST_CHECK_EQUAL(4, std::get<1>(t1));
    std::tuple<Integer, Integer> t2 = ed::divMod(100, -12);
    BOOST_CHECK_EQUAL(-9, std::get<0>(t2));
    BOOST_CHECK_EQUAL(-8, std::get<1>(t2));
    std::tuple<Integer, Integer> t3 = ed::divMod(-100, 12);
    BOOST_CHECK_EQUAL(-9, std::get<0>(t3));
    BOOST_CHECK_EQUAL(8, std::get<1>(t3));
    std::tuple<Integer, Integer> t4 = ed::divMod(-100, -12);
    BOOST_CHECK_EQUAL(8, std::get<0>(t4));
    BOOST_CHECK_EQUAL(-4, std::get<1>(t4));
}

BOOST_AUTO_TEST_CASE(test_rank)
{
    BOOST_CHECK_EQUAL(5_mpz, ed::rank(5_mpz));
    BOOST_CHECK_EQUAL(31268_mpz, ed::rank(ZOmega(1, -2, -9, 12)));
    BOOST_CHECK_EQUAL(4126_mpz, ed::rank(ZRootTwo(-72, 23)));
    BOOST_CHECK_EQUAL(5713_mpz, ed::rank(ZComplex(-72, 23)));
}

BOOST_AUTO_TEST_CASE(test_divmod)
{
    std::tuple<Integer, Integer> t1 = ed::divmod(123_mpz, -24_mpz);
    BOOST_CHECK_EQUAL(-6, std::get<0>(t1));
    BOOST_CHECK_EQUAL(-21, std::get<1>(t1));
    std::tuple<ZRootTwo, ZRootTwo> t2 = ed::divmod(ZRootTwo(-3, -7), ZRootTwo(5, -12));
    BOOST_CHECK_EQUAL(ZRootTwo(1), std::get<0>(t2));
    BOOST_CHECK_EQUAL(ZRootTwo(-8, 5), std::get<1>(t2));
    std::tuple<ZOmega, ZOmega> t3 = ed::divmod(ZOmega(-3, 12, -45, -7), ZOmega(5, -12, 6, 13));
    BOOST_CHECK_EQUAL(ZOmega(-1, 0, -2, -2), std::get<0>(t3));
    BOOST_CHECK_EQUAL(ZOmega(-4, -5, 5, 3), std::get<1>(t3));
    std::tuple<ZComplex, ZComplex> t4 = ed::divmod(ZComplex(-3, -7), ZComplex(5, 12));
    BOOST_CHECK_EQUAL(ZComplex(-1), std::get<0>(t4));
    BOOST_CHECK_EQUAL(ZComplex(2, 5), std::get<1>(t4));
}

BOOST_AUTO_TEST_CASE(test_euclid_mod)
{
    BOOST_CHECK_EQUAL(-21, ed::euclid_mod(123_mpz, Integer(-24)));
    BOOST_CHECK_EQUAL(ZRootTwo(-8, 5), ed::euclid_mod(ZRootTwo(-3, -7), ZRootTwo(5, -12)));
    BOOST_CHECK_EQUAL(ZOmega(-4, -5, 5, 3), ed::euclid_mod(ZOmega(-3, 12, -45, -7), ZOmega(5, -12, 6, 13)));
    BOOST_CHECK_EQUAL(ZComplex(2, 5), ed::euclid_mod(ZComplex(-3, -7), ZComplex(5, 12)));
}

BOOST_AUTO_TEST_CASE(test_euclid_div)
{
    BOOST_CHECK_EQUAL(-6, ed::euclid_div(123_mpz, Integer(-24)));
    BOOST_CHECK_EQUAL(ZRootTwo(1), ed::euclid_div(ZRootTwo(-3, -7), ZRootTwo(5, -12)));
    BOOST_CHECK_EQUAL(ZOmega(-1, 0, -2, -2), ed::euclid_div(ZOmega(-3, 12, -45, -7), ZOmega(5, -12, 6, 13)));
    BOOST_CHECK_EQUAL(ZComplex(-1), ed::euclid_div(ZComplex(-3, -7), ZComplex(5, 12)));
}

BOOST_AUTO_TEST_CASE(test_euclid_gcd)
{
    BOOST_CHECK_EQUAL(-121, ed::euclid_gcd(Integer(26741), Integer(-363)));
    BOOST_CHECK_EQUAL(ZRootTwo(-3, -12), ed::euclid_gcd(ZRootTwo(6, 24), ZRootTwo(-9, -36)));
    BOOST_CHECK_EQUAL(ZOmega(1, 1, 1, 0), ed::euclid_gcd(ZOmega(-3, 12, -45, -7), ZOmega(5, -12, 6, 13)));
    BOOST_CHECK_EQUAL(ZComplex(1), ed::euclid_gcd(ZComplex(-7, 12), ZComplex(12, 12)));
}

BOOST_AUTO_TEST_CASE(test_extended_euclid)
{
    Tup5<Integer> t1 = ed::extended_euclid(Integer(26741), Integer(-363));
    // (1,74,-3,-221,-121)
    BOOST_CHECK_EQUAL(1, std::get<0>(t1));
    BOOST_CHECK_EQUAL(74, std::get<1>(t1));
    BOOST_CHECK_EQUAL(-3, std::get<2>(t1));
    BOOST_CHECK_EQUAL(-221, std::get<3>(t1));
    BOOST_CHECK_EQUAL(-121, std::get<4>(t1));
    Tup5<ZRootTwo> t2 = ed::extended_euclid(ZRootTwo(6, 24), ZRootTwo(-9, -36));
    BOOST_CHECK_EQUAL(ZRootTwo(1), std::get<0>(t2));
    BOOST_CHECK_EQUAL(ZRootTwo(1), std::get<1>(t2));
    BOOST_CHECK_EQUAL(ZRootTwo(-3), std::get<2>(t2));
    BOOST_CHECK_EQUAL(ZRootTwo(-2), std::get<3>(t2));
    BOOST_CHECK_EQUAL(ZRootTwo(-3, -12), std::get<4>(t2));
    Tup5<ZOmega> t3 = ed::extended_euclid(ZOmega(-3, 12, -45, -7), ZOmega(5, -12, 6, 13));
    BOOST_CHECK_EQUAL(ZOmega(-7, 0, 11, -13), std::get<0>(t3));
    BOOST_CHECK_EQUAL(ZOmega(-24, 25, -4, -21), std::get<1>(t3));
    BOOST_CHECK_EQUAL(ZOmega(-5, -12, 30, -23), std::get<2>(t3));
    BOOST_CHECK_EQUAL(ZOmega(-50, 35, 22, -60), std::get<3>(t3));
    BOOST_CHECK_EQUAL(ZOmega(1, 1, 1, 0), std::get<4>(t3));
    Tup5<ZComplex> t4 = ed::extended_euclid(ZComplex(-7, 12), ZComplex(12, 12));
    BOOST_CHECK_EQUAL(ZComplex(5), std::get<0>(t4));
    BOOST_CHECK_EQUAL(ZComplex(-1, -4), std::get<1>(t4));
    BOOST_CHECK_EQUAL(ZComplex(-12, -12), std::get<2>(t4));
    BOOST_CHECK_EQUAL(ZComplex(-7, 12), std::get<3>(t4));
    BOOST_CHECK_EQUAL(ZComplex(1), std::get<4>(t4));
}

BOOST_AUTO_TEST_CASE(test_euclid_inverse)
{
    BOOST_CHECK(std::nullopt == ed::euclid_inverse(Integer(-12)));
    BOOST_CHECK(std::nullopt == ed::euclid_inverse(Integer(12)));
    BOOST_CHECK(1 == ed::euclid_inverse(Integer(1)));
    BOOST_CHECK(-1 == ed::euclid_inverse(Integer(-1)));

    BOOST_CHECK(std::nullopt == ed::euclid_inverse(ZComplex(-12)));
    BOOST_CHECK(std::nullopt == ed::euclid_inverse(ZComplex(12)));
    BOOST_CHECK(ZComplex(0, -1) == ed::euclid_inverse(ZComplex(0, 1)));
    BOOST_CHECK(ZComplex(-1, 0) == ed::euclid_inverse(ZComplex(-1, 0)));

    BOOST_CHECK(std::nullopt == ed::euclid_inverse(ZRootTwo(0, 1)));
    BOOST_CHECK(std::nullopt == ed::euclid_inverse(ZRootTwo(0, -1)));
    BOOST_CHECK(ZRootTwo(1) == ed::euclid_inverse(ZRootTwo(1)));
    BOOST_CHECK(ZRootTwo(-1) == ed::euclid_inverse(ZRootTwo(-1)));

    BOOST_CHECK(std::nullopt == ed::euclid_inverse(ZOmega(0, 0, 2, 0)));
    BOOST_CHECK(std::nullopt == ed::euclid_inverse(ZOmega(1, 0, 1, 0)));
    BOOST_CHECK(ZOmega(-1, 0, 0, 0) == ed::euclid_inverse(ZOmega(0, 0, 1, 0)));
    BOOST_CHECK(ZOmega(0, 1, 0, 0) == ed::euclid_inverse(ZOmega(0, -1, 0, 0)));
}

BOOST_AUTO_TEST_CASE(test_is_unit)
{
    BOOST_CHECK(false == ed::is_unit(Integer(-12)));
    BOOST_CHECK(false == ed::is_unit(Integer(12)));
    BOOST_CHECK(true == ed::is_unit(Integer(1)));
    BOOST_CHECK(true == ed::is_unit(Integer(-1)));
    BOOST_CHECK(false == ed::is_unit(ZComplex(-12)));
    BOOST_CHECK(false == ed::is_unit(ZComplex(12)));
    BOOST_CHECK(true == ed::is_unit(ZComplex(0, 1)));
    BOOST_CHECK(true == ed::is_unit(ZComplex(-1, 0)));
    BOOST_CHECK(false == ed::is_unit(ZRootTwo(0, 1)));
    BOOST_CHECK(false == ed::is_unit(ZRootTwo(0, -1)));
    BOOST_CHECK(true == ed::is_unit(ZRootTwo(1)));
    BOOST_CHECK(true == ed::is_unit(ZRootTwo(-1)));
    BOOST_CHECK(false == ed::is_unit(ZOmega(0, 0, 2, 0)));
    BOOST_CHECK(false == ed::is_unit(ZOmega(1, 0, 1, 0)));
    BOOST_CHECK(true == ed::is_unit(ZOmega(0, 0, 1, 0)));
    BOOST_CHECK(true == ed::is_unit(ZOmega(0, -1, 0, 0)));
}

BOOST_AUTO_TEST_CASE(test_inv_mod)
{
    BOOST_CHECK(std::nullopt == ed::inv_mod(Integer(12), Integer(2)));
    BOOST_CHECK(17 == ed::inv_mod(Integer(22), Integer(13)));
    BOOST_CHECK(std::nullopt == ed::inv_mod(ZRootTwo(12), ZRootTwo(2)));
    BOOST_CHECK(ZRootTwo(-5) == ed::inv_mod(ZRootTwo(22), ZRootTwo(13)));
    BOOST_CHECK(std::nullopt == ed::inv_mod(ZComplex(12), ZComplex(2)));
    BOOST_CHECK(ZComplex(-5) == ed::inv_mod(ZComplex(22), ZComplex(13)));
    BOOST_CHECK(std::nullopt == ed::inv_mod(ZOmega(12), ZOmega(2)));
    BOOST_CHECK(ZOmega(-5) == ed::inv_mod(ZOmega(22), ZOmega(13)));
}

BOOST_AUTO_TEST_CASE(teste_euclid_divides)
{
    BOOST_CHECK(true == ed::euclid_divides(Integer(0), Integer(0)));
    BOOST_CHECK(false == ed::euclid_divides(Integer(0), Integer(9)));
    BOOST_CHECK(true == ed::euclid_divides(Integer(9), Integer(99)));
    BOOST_CHECK(true == ed::euclid_divides(ZRootTwo(0), ZRootTwo(0)));
    BOOST_CHECK(false == ed::euclid_divides(ZRootTwo(0), ZRootTwo(1, 2)));
    BOOST_CHECK(true == ed::euclid_divides(ZRootTwo(1, 2), ZRootTwo(-7, 0)));
    BOOST_CHECK(true == ed::euclid_divides(ZComplex(0), ZComplex(0)));
    BOOST_CHECK(false == ed::euclid_divides(ZComplex(0), ZComplex(1, 2)));
    BOOST_CHECK(true == ed::euclid_divides(ZComplex(1, 2), ZComplex(5, 0)));
    BOOST_CHECK(true == ed::euclid_divides(ZOmega(0), ZOmega(0)));
    BOOST_CHECK(false == ed::euclid_divides(ZOmega(0), ZOmega(1, 2, 3, 4)));
    BOOST_CHECK(true == ed::euclid_divides(ZOmega(1, 2, 3, 4), ZOmega(2, 4, 6, 8)));
}

BOOST_AUTO_TEST_CASE(test_euclid_associates)
{
    BOOST_CHECK(true == ed::euclid_associates(Integer(5), Integer(-5)));
    BOOST_CHECK(false == ed::euclid_associates(Integer(5), Integer(7)));
    BOOST_CHECK(true == ed::euclid_associates(ZRootTwo(5), ZRootTwo(-5)));
    BOOST_CHECK(false == ed::euclid_associates(ZRootTwo(5), ZRootTwo(7)));
    BOOST_CHECK(true == ed::euclid_associates(ZComplex(5), ZComplex(-5)));
    BOOST_CHECK(false == ed::euclid_associates(ZComplex(5), ZComplex(7)));
    BOOST_CHECK(true == ed::euclid_associates(ZOmega(5), ZOmega(-5)));
    BOOST_CHECK(false == ed::euclid_associates(ZOmega(5), ZOmega(7)));
}

BOOST_AUTO_TEST_CASE(test_euclid_extract_power)
{
    BOOST_CHECK(0 == std::get<0>(ed::euclid_extract_power(Integer(0), Integer(12))));
    BOOST_CHECK(0 == std::get<1>(ed::euclid_extract_power(Integer(0), Integer(12))));
    BOOST_CHECK(0 == std::get<0>(ed::euclid_extract_power(Integer(12), Integer(-1))));
    BOOST_CHECK(12 == std::get<1>(ed::euclid_extract_power(Integer(12), Integer(-1))));
    BOOST_CHECK(11 == std::get<0>(ed::euclid_extract_power(Integer(6144), Integer(-2))));
    BOOST_CHECK(-3 == std::get<1>(ed::euclid_extract_power(Integer(6144), Integer(-2))));
    BOOST_CHECK(ZRootTwo(11) == std::get<0>(ed::euclid_extract_power(ZRootTwo(6144), ZRootTwo(-2))));
    BOOST_CHECK(ZRootTwo(-3) == std::get<1>(ed::euclid_extract_power(ZRootTwo(6144), ZRootTwo(-2))));
    BOOST_CHECK(ZComplex(11) == std::get<0>(ed::euclid_extract_power(ZComplex(6144), ZComplex(-2))));
    BOOST_CHECK(ZComplex(-3) == std::get<1>(ed::euclid_extract_power(ZComplex(6144), ZComplex(-2))));
    BOOST_CHECK(ZOmega(11) == std::get<0>(ed::euclid_extract_power(ZOmega(6144), ZOmega(-2))));
    BOOST_CHECK(ZOmega(-3) == std::get<1>(ed::euclid_extract_power(ZOmega(6144), ZOmega(-2))));
}