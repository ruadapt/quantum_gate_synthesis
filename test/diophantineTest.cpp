#define BOOST_TEST_MODULE diophantine
#include "../src/diophantine.h"
#include <boost/test/included/unit_test.hpp>

namespace dio = diophantine;
namespace ed = euclidean_domain;

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
        List<std::tuple<Integer, Integer>> fs;
        std::tie(n, fs) = dio::relatively_prime_factors<Integer>(12, 2);
        BOOST_CHECK_EQUAL(1, n);
        BOOST_CHECK_EQUAL((List<Pair<Integer>>{Pair<Integer>(3, 1), Pair<Integer>(2, 3)}), fs);
    }
    {
        Integer n;
        List<std::tuple<Integer, Integer>> fs;
        std::tie(n, fs) = dio::relatively_prime_factors<Integer>(331500, 6007800);
        BOOST_CHECK_EQUAL(1, n);
        BOOST_CHECK_EQUAL((List<Pair<Integer>>{Pair<Integer>(589, 1), Pair<Integer>(13, 1), Pair<Integer>(5, 5), Pair<Integer>(2, 5), Pair<Integer>(51, 2)}), fs);
    }
    {
        ZOmega z;
        List<std::tuple<ZOmega, Integer>> fs;
        std::tie(z, fs) = dio::relatively_prime_factors<ZOmega>(ZOmega(1, 2, 3, 4), ZOmega(2, 4, 6, 8));
        BOOST_CHECK_EQUAL(ZOmega(0, 0, 0, -1), z);
        std::tuple<ZOmega, Integer> t1 = std::make_tuple(ZOmega(-1, 0, -1, 0), 4);
        std::tuple<ZOmega, Integer> t2 = std::make_tuple(ZOmega(3, 1, 1, -2), 2);
        BOOST_CHECK((List<std::tuple<ZOmega, Integer>>{t1, t2}) == fs);
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
    Integer val = s.run();
    BOOST_CHECK(8 == val || 5 == val);
}

BOOST_AUTO_TEST_CASE(test_root_mod)
{
    {
        StepComp<Integer> s = dio::root_mod(13, 4);
        Integer val = s.run();
        BOOST_CHECK(2 == val || 11 == val);
    }
    {
        StepComp<Integer> s = dio::root_mod(129, 13);
        Integer val = s.run();
        BOOST_CHECK(20 == val || 23 == val || 106 == val || 109 == val);
    }
}

BOOST_AUTO_TEST_CASE(test_dioph_int_assoc_prime)
{
    // n < 0 case
    {
        StepComp<Maybe<ZOmega>> sc = dio::dioph_int_assoc_prime(-11);
        Maybe<ZOmega> z = sc.run();
        BOOST_CHECK(ZOmega(1, 0, 1, 3) == z || ZOmega(1, 0, 1, -3) == z);
    }
    // n == 0 case
    {
        StepComp<Maybe<ZOmega>> sc = dio::dioph_int_assoc_prime(0);
        Maybe<ZOmega> z = sc.run();
        BOOST_CHECK(ZOmega(0) == z);
    }
    // n == 2 case
    {
        StepComp<Maybe<ZOmega>> sc = dio::dioph_int_assoc_prime(2);
        Maybe<ZOmega> z = sc.run();
        BOOST_CHECK(ZOmega(-1, 0, 1, 0) == z);
    }
    // n mod 4 == 1 case
    {
        StepComp<Maybe<ZOmega>> sc = dio::dioph_int_assoc_prime(5);
        Maybe<ZOmega> z = sc.run();
        BOOST_CHECK(ZOmega(0, 1, 0, 2) == z || ZOmega(0, 1, 0, -2) == z);
    }
    // n mod 8 == 3 case
    {
        StepComp<Maybe<ZOmega>> sc = dio::dioph_int_assoc_prime(11);
        Maybe<ZOmega> z = sc.run();
        BOOST_CHECK(ZOmega(1, 0, 1, 3) == z || ZOmega(1, 0, 1, -3) == z);
    }
    // n mod 8 == 7 case
    {
        StepComp<Maybe<ZOmega>> sc = dio::dioph_int_assoc_prime(23);
        Maybe<ZOmega> z = sc.run();
        BOOST_CHECK(Maybe<ZOmega>() == z);
    }
}

BOOST_AUTO_TEST_CASE(test_dioph_int_assoc)
{
    {
        StepComp<Maybe<ZOmega>> sc = dio::dioph_int_assoc(171);
        Maybe<ZOmega> mz = sc.run();
        BOOST_REQUIRE(mz.has_value());
        ZOmega z = mz.value();
        ZOmega prod = z.adj() * z;
        // Make sure that the returned value is a correct solution.
        BOOST_CHECK(ed::euclid_associates(ZOmega(171), prod));
    }
    {
        StepComp<Maybe<ZOmega>> sc = dio::dioph_int_assoc(1235);
        Maybe<ZOmega> mz = sc.run();
        BOOST_REQUIRE(mz.has_value());
        ZOmega z = mz.value();
        ZOmega prod = z.adj() * z;
        // Make sure that the returned value is a correct solution.
        BOOST_CHECK(ed::euclid_associates(ZOmega(1235), prod));
    }
    {
        StepComp<Maybe<ZOmega>> sc = dio::dioph_int_assoc(12367849);
        Maybe<ZOmega> mz = sc.run();
        BOOST_REQUIRE(mz.has_value());
        ZOmega z = mz.value();
        ZOmega prod = z.adj() * z;
        // Make sure that the returned value is a correct solution.
        BOOST_CHECK(ed::euclid_associates(ZOmega(12367849), prod));
    }
}

BOOST_AUTO_TEST_CASE(test_diophantine)
{
    // has_result indicates whether the diophantine equation should have a solution.
    auto test = [](Integer a, Integer b, bool has_result)
    {
        StepComp<Maybe<ZOmega>> sc = dio::diophantine(ZRootTwo(a, b));
        Maybe<ZOmega> mz = sc.run();
        if (!has_result)
        {
            BOOST_REQUIRE(!mz.has_value());
            return;
        }
        BOOST_REQUIRE(mz.has_value());
        ZOmega z = mz.value();
        ZOmega prod = z.adj() * z;
        ZOmega expected_prod = ZOmega(a) + ZOmega(b) * ring::roottwo<ZOmega>();
        BOOST_CHECK_EQUAL(expected_prod, prod);
    };
    test(6_mpz, 1_mpz, true);
    test(18_mpz, -11_mpz, true);
    test(1230948_mpz, 873216_mpz, false);
}
