#include "quadratic.h"
#include "toReal.h"
#include <cmath>

template <typename T>
std::optional<std::tuple<Real, Real>> quadratic(T a_, T b_, T c_)
{
    Real a = toReal<Real>(a_);
    Real b = toReal<Real>(b_);
    Real c = toReal<Real>(c_);
    Real radix = b * b - 4 * a * c;
    if (radix < 0)
    {
        return std::nullopt;
    }
    Real s1 = -b - sqrt(radix);
    Real s2 = -b + sqrt(radix);
    if (b >= 0)
    {
        Real t1 = s1 / (2 * a);
        Real t2 = (2 * c) / s1;
        return std::make_tuple(t1, t2);
    }
    Real t1prime = (2 * c) / s2;
    Real t2prime = s2 / (2 * a);
    return std::make_tuple(t1prime, t2prime);
}