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