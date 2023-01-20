#define BOOST_TEST_MODULE toReal
#include "comparisons.h"
#include "../toReal.h"
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE(test_simple_types)
{
    BOOST_TEST(1.25 == toReal<Real>(5_mpq / 4));
    BOOST_TEST(12.0 == toReal<Real>(12_mpz));
    BOOST_TEST(12.0 == toReal<Real>(12));
    // BOOST_TEST(1.56 == toReal<Real>(1.56)); // TODO look back at this
}

BOOST_AUTO_TEST_CASE(test_root_two)
{
    BOOST_TEST(approx_equal(2.414213562373095, toReal<Real>(ZRootTwo(1, 1))));
    BOOST_TEST(approx_equal(3.4393398282201786, toReal<Real>(DRootTwo(ZDyadic(9, 1), ZDyadic(-3, 2)))));
    BOOST_TEST(approx_equal(1.3285911063657792, toReal<Real>(QRootTwo(12_mpq / 7, -3_mpq / 11))));
}