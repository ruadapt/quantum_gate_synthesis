#include "quadratic.h"
#include "types.h"
#include "toReal.h"
#include "ring.h"

template <typename T>
Maybe<std::tuple<Real, Real>> quadratic(T a_, T b_, T c_)
{
    Real a = to_real<Real>(a_);
    Real b = to_real<Real>(b_);
    Real c = to_real<Real>(c_);
    // This check isn't in the Haskell version, but without it we get NaN
    // in this case.
    if (b == 0 && c == 0)
    {
        return std::make_tuple(0, 0);
    }
    Real radix = b * b - 4 * a * c;
    if (radix < 0)
    {
        return std::nullopt;
    }
    Real s1 = -b - bmp::sqrt(radix);
    Real s2 = -b + bmp::sqrt(radix);
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