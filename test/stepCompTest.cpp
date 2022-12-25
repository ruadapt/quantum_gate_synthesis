#define BOOST_TEST_MODULE stepComp
#include "../stepComp.h"
#include "../types.h"
#include <boost/test/included/unit_test.hpp>

namespace sc = stepcomp;

BOOST_AUTO_TEST_CASE(test_construction)
{
    StepComp<Integer> s1 = StepComp<Integer>(Integer(2));
    BOOST_CHECK_EQUAL(true, s1.is_done());
    BOOST_CHECK_EQUAL(2, s1.value());

    StepComp<Integer> s2 = StepComp<Integer>([]()
                                             { return StepComp<Integer>(Integer(5)); });
    BOOST_CHECK_EQUAL(false, s2.is_done());
    BOOST_CHECK_EQUAL(5, s2.comp()().value());
}

BOOST_AUTO_TEST_CASE(test_untick)
{
    StepComp<Integer> s1 = StepComp<Integer>(Integer(2));
    BOOST_CHECK_EQUAL(2, s1.untick().value());

    StepComp<Integer> s2 = StepComp<Integer>([]()
                                             { return StepComp<Integer>(Integer(5)); });
    BOOST_CHECK_EQUAL(false, s2.is_done());
    BOOST_CHECK_EQUAL(true, s2.untick().is_done());
    BOOST_CHECK_EQUAL(5, s2.untick().value());

    StepComp<Integer> s3 = sc::wrap(Integer(7), 2);
    BOOST_CHECK_EQUAL(false, s3.is_done());
    BOOST_CHECK_EQUAL(false, s3.untick().is_done());
    BOOST_CHECK_EQUAL(true, s3.untick().untick().is_done());
    BOOST_CHECK_EQUAL(7, s3.untick().untick().value());
}

BOOST_AUTO_TEST_CASE(test_forward)
{
    StepComp<Integer> s1 = StepComp(Integer(10));
    BOOST_CHECK_EQUAL(true, s1.forward(100).is_done());
    BOOST_CHECK_EQUAL(10, s1.forward(100).value());

    StepComp<Integer> s2 = sc::wrap(Integer(15), 100);
    BOOST_CHECK_EQUAL(false, s2.is_done());
    BOOST_CHECK_EQUAL(false, s2.forward(99).is_done());
    BOOST_CHECK_EQUAL(true, s2.forward(100).is_done());
    BOOST_CHECK_EQUAL(15, s2.forward(100).value());
}

BOOST_AUTO_TEST_CASE(test_speedup)
{
    StepComp<Integer> s = sc::wrap(Integer(15), 8);
    BOOST_CHECK_EQUAL(1, s.speed());
    BOOST_CHECK_EQUAL(false, s.is_done());
    BOOST_CHECK_EQUAL(false, s.forward(7).is_done());
    BOOST_CHECK_EQUAL(true, s.forward(8).is_done());

    StepComp<Integer> s2 = s.speedup(2);
    BOOST_CHECK_EQUAL(2, s2.speed());
    BOOST_CHECK_EQUAL(false, s2.is_done());
    BOOST_CHECK_EQUAL(false, s2.forward(3).is_done());
    BOOST_CHECK_EQUAL(true, s2.forward(4).is_done());

    StepComp<Integer> s3 = s2.speedup(2);
    BOOST_CHECK_EQUAL(4, s3.speed());
    BOOST_CHECK_EQUAL(false, s3.is_done());
    BOOST_CHECK_EQUAL(false, s3.forward(1).is_done());
    BOOST_CHECK_EQUAL(true, s3.forward(2).is_done());
}

BOOST_AUTO_TEST_CASE(test_parallel)
{
    StepComp<std::string> sc1 = sc::wrap(std::string("result1"), 3);
    StepComp<std::string> sc2 = sc::wrap(std::string("result2"), 2);
    using A = std::string;
    using B = std::string;
    using R = Either<std::tuple<A, StepComp<B>>, std::tuple<StepComp<A>, B>>;
    StepComp<R> p = sc::parallel(sc1, sc2);
    BOOST_CHECK_EQUAL(false, p.is_done());
    BOOST_CHECK_EQUAL(false, p.forward(1).is_done());
    BOOST_CHECK_EQUAL(true, p.forward(2).is_done());
    StepComp<R> d = p.forward(2);
    BOOST_CHECK_EQUAL(2, d.value().index()); // The second computation finishes first.
    std::string result = std::get<1>(snd(d.value()));
    BOOST_CHECK_EQUAL("result2", result);
}

BOOST_AUTO_TEST_CASE(test_parallel_first)
{
    StepComp<std::string> sc1 = sc::wrap(std::string("result1"), 2);
    StepComp<std::string> sc2 = sc::wrap(std::string("result2"), 3);
    StepComp<std::string> p = sc::parallel_first(sc1, sc2);
    BOOST_CHECK_EQUAL(false, p.is_done());
    BOOST_CHECK_EQUAL(false, p.forward(1).is_done());
    BOOST_CHECK_EQUAL(true, p.forward(2).is_done());
    BOOST_CHECK_EQUAL("result1", p.forward(2).value()); // The first computation finishes first.
}

BOOST_AUTO_TEST_CASE(test_parallel_maybe)
{
    {
        StepComp<Maybe<std::string>> sc1 = sc::wrap<Maybe<std::string>>(std::nullopt, 5);
        StepComp<Maybe<int>> sc2 = sc::wrap<Maybe<int>>(std::nullopt, 4);
        StepComp<Maybe<std::tuple<std::string, int>>> p = sc::parallel_maybe(sc1, sc2);
        BOOST_CHECK_EQUAL(false, p.is_done());
        BOOST_CHECK_EQUAL(false, p.forward(3).is_done());
        BOOST_CHECK_EQUAL(true, p.forward(4).is_done());
        BOOST_CHECK(std::nullopt == p.run());
    }
    {
        StepComp<Maybe<std::string>> sc1 = sc::wrap<Maybe<std::string>>("r1", 5);
        StepComp<Maybe<int>> sc2 = sc::wrap<Maybe<int>>(std::nullopt, 4);
        StepComp<Maybe<std::tuple<std::string, int>>> p = sc::parallel_maybe(sc1, sc2);
        BOOST_CHECK_EQUAL(false, p.is_done());
        BOOST_CHECK_EQUAL(false, p.forward(3).is_done());
        BOOST_CHECK_EQUAL(true, p.forward(4).is_done());
        BOOST_CHECK(std::nullopt == p.run());
    }
    {
        StepComp<Maybe<std::string>> sc1 = sc::wrap<Maybe<std::string>>(std::nullopt, 5);
        StepComp<Maybe<int>> sc2 = sc::wrap<Maybe<int>>(12, 4);
        StepComp<Maybe<std::tuple<std::string, int>>> p = sc::parallel_maybe(sc1, sc2);
        BOOST_CHECK_EQUAL(false, p.is_done());
        BOOST_CHECK_EQUAL(false, p.forward(4).is_done());
        BOOST_CHECK_EQUAL(true, p.forward(5).is_done());
        BOOST_CHECK(std::nullopt == p.run());
    }
    {
        StepComp<Maybe<std::string>> sc1 = sc::wrap<Maybe<std::string>>("r1", 5);
        StepComp<Maybe<int>> sc2 = sc::wrap<Maybe<int>>(12, 4);
        StepComp<Maybe<std::tuple<std::string, int>>> p = sc::parallel_maybe(sc1, sc2);
        BOOST_CHECK_EQUAL(false, p.is_done());
        BOOST_CHECK_EQUAL(false, p.forward(4).is_done());
        BOOST_CHECK_EQUAL(true, p.forward(5).is_done());
        BOOST_CHECK(std::make_tuple("r1", 12) == p.run());
    }
}

BOOST_AUTO_TEST_CASE(test_diverge)
{
    StepComp<Integer> div = sc::diverge<Integer>();
    BOOST_CHECK_EQUAL(false, div.is_done());
    BOOST_CHECK_EQUAL(false, div.forward(1000).is_done());
}

BOOST_AUTO_TEST_CASE(test_count)
{
    StepComp<int> sc = sc::wrap(5, 100);
    StepComp<int> sc2 = sc.speedup(50);

    StepComp<int> f1 = sc.forward(99);
    BOOST_CHECK_EQUAL(false, f1.is_done());

    StepComp<int> f2 = sc.forward(100);
    BOOST_CHECK_EQUAL(true, f2.is_done());
    BOOST_CHECK_EQUAL(100, f2.count());

    // Once the StepComp is done after 100 steps, we stop counting.
    StepComp<int> f3 = sc.forward(300);
    BOOST_CHECK_EQUAL(true, f3.is_done());
    BOOST_CHECK_EQUAL(100, f3.count());

    StepComp<int> f4 = sc2.forward(1);
    BOOST_CHECK_EQUAL(false, f4.is_done());

    // Even with a speedup, we still count the raw number of computations.
    StepComp<int> f5 = sc2.forward(2);
    BOOST_CHECK_EQUAL(true, f5.is_done());
    BOOST_CHECK_EQUAL(100, f5.count());
}