#define BOOST_TEST_MODULE quadratic
#include "comparisons.h"
#include "../src/types.h"
#include "../src/quadratic.h"
#include <boost/test/included/unit_test.hpp>
#include <optional>
#include <tuple>

BOOST_AUTO_TEST_CASE(test_simple_case)
{
    Integer a = 1;
    Integer b = 2;
    Integer c = 1;
    std::optional<std::tuple<Real, Real>> rOpt = quadratic<Integer>(a, b, c);
    BOOST_REQUIRE(rOpt.has_value());
    std::tuple<Real, Real> r = rOpt.value();
    Real s1 = std::get<0>(r);
    Real s2 = std::get<1>(r);
    BOOST_CHECK_EQUAL(-1, s1);
    BOOST_CHECK_EQUAL(-1, s2);
}

BOOST_AUTO_TEST_CASE(test_no_solution)
{
    Integer a = 7;
    Integer b = 0;
    Integer c = 1;
    std::optional<std::tuple<Real, Real>> rOpt = quadratic<Integer>(a, b, c);
    BOOST_REQUIRE(!rOpt.has_value());
}

BOOST_AUTO_TEST_CASE(test_non_integral1)
{
    Integer a = 5;
    Integer b = 13;
    Integer c = -9;
    std::optional<std::tuple<Real, Real>> rOpt = quadratic<Integer>(a, b, c);
    BOOST_REQUIRE(rOpt.has_value());
    std::tuple<Real, Real> r = rOpt.value();
    Real s1 = std::get<0>(r);
    Real s2 = std::get<1>(r);
    Real e1 = -3.1681541692269404;
    Real e2 = 0.5681541692269404;
    BOOST_TEST(approx_equal(e1, s1));
    BOOST_TEST(approx_equal(e2, s2));
}

BOOST_AUTO_TEST_CASE(test_non_integral2)
{
    Integer a = 5;
    Integer b = -13;
    Integer c = -9;
    std::optional<std::tuple<Real, Real>> rOpt = quadratic<Integer>(a, b, c);
    BOOST_REQUIRE(rOpt.has_value());
    std::tuple<Real, Real> r = rOpt.value();
    Real s1 = std::get<0>(r);
    Real s2 = std::get<1>(r);
    Real e1 = -0.5681541692269404;
    Real e2 = 3.1681541692269404;
    BOOST_TEST(approx_equal(e1, s1));
    BOOST_TEST(approx_equal(e2, s2));
}

BOOST_AUTO_TEST_CASE(test_b_and_c_are_zero_case)
{
    QRootTwo a = QRootTwo(3, 2);
    QRootTwo b = 0;
    QRootTwo c = 0;
    Maybe<Pair<Real>> q = quadratic<QRootTwo>(a, b, c);
    BOOST_REQUIRE(q.has_value());
    Pair<Real> expected = {0, 0};
    BOOST_CHECK(expected == q.value());
}

BOOST_AUTO_TEST_CASE(test_b_is_zero_case)
{
    QRootTwo a = 1;
    QRootTwo b = 0;
    QRootTwo c = -4;
    Maybe<Pair<Real>> q = quadratic<QRootTwo>(a, b, c);
    BOOST_REQUIRE(q.has_value());
    Pair<Real> expected = {-2, 2};
    BOOST_CHECK(expected == q.value());
}

BOOST_AUTO_TEST_CASE(test_c_is_zero_case)
{
    QRootTwo a = 1;
    QRootTwo b = 1;
    QRootTwo c = 0;
    Maybe<Pair<Real>> q = quadratic<QRootTwo>(a, b, c);
    BOOST_REQUIRE(q.has_value());
    Pair<Real> expected = {-1, 0};
    BOOST_CHECK(expected == q.value());
}