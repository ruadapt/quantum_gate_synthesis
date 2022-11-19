#include "toReal.h"
#include "ring.h"

template <typename T>
T toReal(Rational r)
{
    return ring::fromRational<T>(r);
}

template <typename T>
T toReal(Integer n)
{
    return ring::fromInteger<T>(n);
}

template <typename T>
T toReal(int n)
{
    return toReal<T>(Integer(n));
}

template <typename T>
T toReal(double d)
{
    return ring::fromRational<T>(Rational(d));
}