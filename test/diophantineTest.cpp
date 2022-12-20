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
        std::cout << "factor = " << factor << std::endl;
        BOOST_CHECK(7907 == factor || 7919 == factor);
    }
}