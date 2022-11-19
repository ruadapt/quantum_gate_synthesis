#include "quadratic.h"
#include "toReal.h"
#include <cmath>

template <typename T>
std::optional<std::tuple<double, double>> quadratic(T a_, T b_, T c_)
{
    double a = toReal<double>(a_);
    double b = toReal<double>(b_);
    double c = toReal<double>(c_);
    double radix = b * b - 4 * a * c;
    if (radix < 0)
    {
        return std::nullopt;
    }
    double s1 = -b - sqrt(radix);
    double s2 = -b + sqrt(radix);
    if (b >= 0)
    {
        double t1 = s1 / (2 * a);
        double t2 = (2 * c) / s1;
        return std::make_tuple(t1, t2);
    }
    double t1prime = (2 * c) / s2;
    double t2prime = s2 / (2 * a);
    return std::make_tuple(t1prime, t2prime);
}