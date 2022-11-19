#define BOOST_TEST_MODULE quadratic
#include "../quadratic.h"
#include <boost/test/included/unit_test.hpp>
#include <optional>
#include <tuple>

BOOST_AUTO_TEST_CASE(simple_case)
{
    Integer a = 1;
    Integer b = 2;
    Integer c = 1;
    std::optional<std::tuple<double, double>> rOpt = quadratic<Integer>(a, b, c);
    BOOST_REQUIRE(rOpt.has_value());
    std::tuple<double, double> r = rOpt.value();
    double s1 = std::get<0>(r);
    double s2 = std::get<1>(r);
    BOOST_CHECK_EQUAL(-1, s1);
    BOOST_CHECK_EQUAL(-1, s2);
}

BOOST_AUTO_TEST_CASE(no_solution)
{
    Integer a = 7;
    Integer b = 0;
    Integer c = 1;
    std::optional<std::tuple<double, double>> rOpt = quadratic<Integer>(a, b, c);
    BOOST_REQUIRE(!rOpt.has_value());
}

BOOST_AUTO_TEST_CASE(non_integral1)
{
    Integer a = 5;
    Integer b = 13;
    Integer c = -9;
    std::optional<std::tuple<double, double>> rOpt = quadratic<Integer>(a, b, c);
    BOOST_REQUIRE(rOpt.has_value());
    std::tuple<double, double> r = rOpt.value();
    double s1 = std::get<0>(r);
    double s2 = std::get<1>(r);
    double e1 = -3.1681541692269404;
    double e2 = 0.5681541692269404;
    BOOST_TEST(e1 == s1, boost::test_tools::tolerance(0.001));
    BOOST_TEST(e2 == s2, boost::test_tools::tolerance(0.001));
}

BOOST_AUTO_TEST_CASE(non_integral2)
{
    Integer a = 5;
    Integer b = -13;
    Integer c = -9;
    std::optional<std::tuple<double, double>> rOpt = quadratic<Integer>(a, b, c);
    BOOST_REQUIRE(rOpt.has_value());
    std::tuple<double, double> r = rOpt.value();
    double s1 = std::get<0>(r);
    double s2 = std::get<1>(r);
    double e1 = -0.5681541692269404;
    double e2 = 3.1681541692269404;
    BOOST_TEST(e1 == s1, boost::test_tools::tolerance(0.001));
    BOOST_TEST(e2 == s2, boost::test_tools::tolerance(0.001));
}