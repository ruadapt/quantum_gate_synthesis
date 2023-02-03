/** \file euclideanDomain.h
 */
#pragma once
#include "ring.h"

namespace euclidean_domain
{
    Integer rank(Integer x);

    Integer rank(ZOmega x);

    Integer rank(ZComplex x);

    Integer rank(ZRootTwo x);

    std::tuple<Integer, Integer> divmod(Integer x, Integer y);

    std::tuple<ZOmega, ZOmega> divmod(ZOmega x, ZOmega y);

    std::tuple<ZComplex, ZComplex> divmod(ZComplex x, ZComplex y);

    std::tuple<ZRootTwo, ZRootTwo> divmod(ZRootTwo x, ZRootTwo y);

    template <typename T>
    T euclid_mod(T x, T y);
    
    template <typename T>
    T euclid_div(T x, T y);

    template <typename T>
    T euclid_gcd(T x, T y);

    template <typename T>
    std::tuple<T, T, T, T, T> extended_euclid(T x, T y);

    template <typename T>
    Maybe<T> euclid_inverse(T x);

    template <typename T>
    bool is_unit(T x);

    template <typename T>
    Maybe<T> inv_mod(T p, T a);

    template <typename T>
    bool euclid_divides(T a, T b);

    template <typename T>
    bool euclid_associates(T a, T b);

    template <typename T>
    std::tuple<Integer, T> euclid_extract_power(T x, T y);

    Integer rounddiv(Integer x, Integer y);

    std::tuple<Integer, Integer> divMod(Integer x, Integer y);
}

#include "euclideanDomain.cpp"