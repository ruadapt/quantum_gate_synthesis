#define BOOST_TEST_MODULE toReal
#include "comparisons.h"
#include "../src/toReal.h"
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE(test_simple_types)
{
    BOOST_TEST(1.25 == to_real<Real>(5_mpq / 4));
    BOOST_TEST(12.0 == to_real<Real>(12_mpz));
    BOOST_TEST(12.0 == to_real<Real>(12));
}

BOOST_AUTO_TEST_CASE(test_root_two)
{
    BOOST_TEST(approx_equal(2.414213562373095, to_real<Real>(ZRootTwo(1, 1))));
    BOOST_TEST(approx_equal(3.4393398282201786, to_real<Real>(DRootTwo(ZDyadic(9, 1), ZDyadic(-3, 2)))));
    BOOST_TEST(approx_equal(1.3285911063657792, to_real<Real>(QRootTwo(12_mpq / 7, -3_mpq / 11))));
}