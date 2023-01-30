/** \file toReal.cpp
 */
#include "toReal.h"
#include "ring.h"

template <typename T>
T to_real(Rational r)
{
    return ring::fromRational<T>(r);
}

template <typename T>
T to_real(Integer n)
{
    return ring::fromInteger<T>(n);
}

template <typename T>
T to_real(int n)
{
    return to_real<T>(Integer(n));
}

// TODO see if this is needed
// template <typename T>
// T to_real(Real r)
// {
//     return ring::fromRational<T>(Rational(r));
// }

template <typename T>
T to_real(ZDyadic d)
{
    return to_real<T>(d.a()) / to_real<T>(ring::exp2<Integer>(d.n()));
}

template <typename T>
T to_real(ZRootTwo r)
{
    return to_real<T>(r.a()) + ring::roottwo<T>() * to_real<T>(r.b());
}

template <typename T>
T to_real(DRootTwo r)
{
    return to_real<T>(r.a()) + ring::roottwo<T>() * to_real<T>(r.b());
}

template <typename T>
T to_real(QRootTwo r)
{
    return to_real<T>(r.a()) + ring::roottwo<T>() * to_real<T>(r.b());
}