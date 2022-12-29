#define BOOST_TEST_MODULE diophantine
#include "../diophantine.h"
#include <boost/test/included/unit_test.hpp>

namespace dio = diophantine;

BOOST_AUTO_TEST_CASE(test_find_factor)
{
    {
        StepComp<Integer> s = dio::find_factor(15);
        Integer factor = s.run();
        BOOST_CHECK(3 == factor || 5 == factor);
    }
    {
        StepComp<Integer> s = dio::find_factor(62615533); // 7907 * 7919
        Integer factor = s.run();
        BOOST_CHECK(7907 == factor || 7919 == factor);
    }
    {
        StepComp<Integer> s = dio::find_factor(373782744163); // 610817 * 611939
        Integer factor = s.run();
        BOOST_CHECK(610817 == factor || 611939 == factor);
    }
}

BOOST_AUTO_TEST_CASE(test_relatively_prime_factors)
{
    {
        Integer n;
        List<Pair<Integer>> fs;
        std::tie(n, fs) = dio::relatively_prime_factors(12, 2);
        BOOST_CHECK_EQUAL(1, n);
        BOOST_CHECK_EQUAL((List<Pair<Integer>>{Pair<Integer>(3, 1), Pair<Integer>(2, 3)}), fs);
    }
    {
        Integer n;
        List<Pair<Integer>> fs;
        std::tie(n, fs) = dio::relatively_prime_factors(331500, 6007800);
        BOOST_CHECK_EQUAL(1, n);
        BOOST_CHECK_EQUAL((List<Pair<Integer>>{Pair<Integer>(589, 1), Pair<Integer>(13, 1), Pair<Integer>(5, 5), Pair<Integer>(2, 5), Pair<Integer>(51, 2)}), fs);
    }
}

BOOST_AUTO_TEST_CASE(test_power_mod)
{
    BOOST_CHECK_EQUAL(1, dio::power_mod(2, 0, 12));
    BOOST_CHECK_EQUAL(1, dio::power_mod(234, 0, 17));
    BOOST_CHECK_EQUAL(5, dio::power_mod(17, 1, 12));
    BOOST_CHECK_EQUAL(6, dio::power_mod(17, 123, 7));
}

BOOST_AUTO_TEST_CASE(test_root_of_negative_one)
{
    StepComp<Integer> s = dio::root_of_negative_one(13);
    BOOST_CHECK(8 == s.run() || 5 == s.run());
}