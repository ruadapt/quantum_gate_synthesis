#define BOOST_TEST_MODULE toReal
#include "../toReal.h"
#include <boost/test/included/unit_test.hpp>

namespace tt = boost::test_tools;

BOOST_AUTO_TEST_CASE(simple_types)
{
    BOOST_TEST(1.25 == toReal<Real>(5_mpq / 4));
    BOOST_TEST(12.0 == toReal<Real>(12_mpz));
    BOOST_TEST(12.0 == toReal<Real>(12));
    BOOST_TEST(1.56 == toReal<Real>(1.56));
}

BOOST_AUTO_TEST_CASE(root_two)
{
    BOOST_TEST(2.414 == toReal<Real>(ZRootTwo(1, 1)), tt::tolerance(0.001));
    BOOST_TEST(3.439 == toReal<Real>(DRootTwo(ZDyadic(9, 1), ZDyadic(-3, 2))), tt::tolerance(0.001));
    BOOST_TEST(1.3286 == toReal<Real>(QRootTwo(12_mpq / 7, -3_mpq / 11)), tt::tolerance(0.001));
}