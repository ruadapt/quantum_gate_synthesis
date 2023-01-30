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

// TODO see if this is needed
// template <typename T>
// T toReal(Real r)
// {
//     return ring::fromRational<T>(Rational(r));
// }

template <typename T>
T toReal(ZDyadic d)
{
    return toReal<T>(d.a()) / toReal<T>(ring::exp2<Integer>(d.n()));
}

template <typename T>
T toReal(ZRootTwo r)
{
    return toReal<T>(r.a()) + ring::rootTwo<T>() * toReal<T>(r.b());
}

template <typename T>
T toReal(DRootTwo r)
{
    return toReal<T>(r.a()) + ring::rootTwo<T>() * toReal<T>(r.b());
}

template <typename T>
T toReal(QRootTwo r)
{
    return toReal<T>(r.a()) + ring::rootTwo<T>() * toReal<T>(r.b());
}