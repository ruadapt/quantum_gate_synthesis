#include "../src/types.h"

bool approx_equal(Real x1, Real x2)
{
    return boost::multiprecision::abs(x1 - x2) < 0.0001;
}

bool approx_equal(double x1, double x2)
{
    return abs(x1 - x2) < 0.0001;
}