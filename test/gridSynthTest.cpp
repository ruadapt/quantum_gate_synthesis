#define BOOST_TEST_MODULE gridSynth
#include "comparisons.h"
#include "../src/gridSynth.h"
#include <boost/test/included/unit_test.hpp>

namespace gs = gridsynth;

BOOST_AUTO_TEST_CASE(test_epsilon_region_creation)
{
    ConvexSet<Real> eps = gs::epsilon_region<Real>(0.2, 0.4);
}

BOOST_AUTO_TEST_CASE(test_gridsynth_internal)
{
    Real prec = 5;
    Real theta = 1;
    int effort = 25;

    U2<DOmega> uU;
    Maybe<double> error;
    List<std::tuple<DOmega, Integer, gs::DStatus>> candidate_info;
    std::tie(uU, error, candidate_info) = gs::gridsynth_internal<Real>(prec, theta, effort);

    DOmega expected_a = DOmega(ZDyadic(-9, 4), ZDyadic(-7, 4), ZDyadic(1, 1), ZDyadic(1, 3));
    DOmega expected_b = DOmega(ZDyadic(-1, 3), ZDyadic(5, 4), ZDyadic(-5, 4), ZDyadic(1, 3));
    DOmega expected_c = DOmega(ZDyadic(-5, 4), ZDyadic(5, 4), ZDyadic(-1, 3), ZDyadic(-1, 3));
    DOmega expected_d = DOmega(ZDyadic(-1, 1), ZDyadic(7, 4), ZDyadic(9, 4), ZDyadic(1, 3));

    BOOST_CHECK_EQUAL(expected_a, uU(0, 0));
    BOOST_CHECK_EQUAL(expected_b, uU(0, 1));
    BOOST_CHECK_EQUAL(expected_c, uU(1, 0));
    BOOST_CHECK_EQUAL(expected_d, uU(1, 1));
}

BOOST_AUTO_TEST_CASE(test_gridsynth_internal2)
{
    Real prec = 5;
    Real theta = 1.5707963267948966;
    int effort = 25;

    U2<DOmega> uU;
    Maybe<double> error;
    List<std::tuple<DOmega, Integer, gs::DStatus>> candidate_info;
    std::tie(uU, error, candidate_info) = gs::gridsynth_internal<Real>(prec, theta, effort);
    
    DOmega expected_a = DOmega(-1, 0, 0, 0);
    DOmega expected_b = DOmega(0, 0, 0, 0);
    DOmega expected_c = DOmega(0, 0, 0, 0);
    DOmega expected_d = DOmega(0, 0, 1, 0);

    BOOST_CHECK_EQUAL(expected_a, uU(0, 0));
    BOOST_CHECK_EQUAL(expected_b, uU(0, 1));
    BOOST_CHECK_EQUAL(expected_c, uU(1, 0));
    BOOST_CHECK_EQUAL(expected_d, uU(1, 1));
}

BOOST_AUTO_TEST_CASE(test_gridsynth_internal3)
{
    Real prec = 5;
    Real theta = 1.2;
    int effort = 25;

    U2<DOmega> uU;
    Maybe<double> error;
    List<std::tuple<DOmega, Integer, gs::DStatus>> candidate_info;
    std::tie(uU, error, candidate_info) = gs::gridsynth_internal<Real>(prec, theta, effort);

    DOmega expected_a = DOmega(ZDyadic(-1, 1), ZDyadic(-5, 3), ZDyadic(9, 4), ZDyadic(1, 4));
    DOmega expected_b = DOmega(ZDyadic(1, 4), ZDyadic(-1, 3), ZDyadic(1, 3), ZDyadic(-1, 4));
    DOmega expected_c = DOmega(ZDyadic(1, 3), ZDyadic(-1, 3), ZDyadic(1, 4), ZDyadic(1, 4));
    DOmega expected_d = DOmega(ZDyadic(-9, 4), ZDyadic(5, 3), ZDyadic(1, 1), ZDyadic(1, 4));

    BOOST_CHECK_EQUAL(expected_a, uU(0, 0));
    BOOST_CHECK_EQUAL(expected_b, uU(0, 1));
    BOOST_CHECK_EQUAL(expected_c, uU(1, 0));
    BOOST_CHECK_EQUAL(expected_d, uU(1, 1));
}