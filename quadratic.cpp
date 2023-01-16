#include "quadratic.h"
#include "types.h"
#include "toReal.h"
#include "ring.h"

template <typename T>
Maybe<Pair<Integer>> int_quadratic(T b, T c)
{
    T radix = (b * b) * ring::recip<T>(T(4)) - c;
    if (radix < 0)
    {
        return Maybe<Pair<Integer>>();
    }
    T tm = -b * ring::recip<T>(T(2));
    Integer rootradix_prime = ring::intsqrt(ring::floor_of(radix));
    auto f = [b, c](T x) -> T
    {
        return x * x + b * x + c;
    };
    auto is_solution1 = [f, tm](Integer x) -> bool
    {
        T x_prime = ring::fromInteger<T>(x);
        return f(x_prime) >= 0 && (f(x_prime - 1) < 0 || x_prime - 1 < tm);
    };
    auto is_solution0 = [f, tm](Integer x) -> bool
    {
        T x_prime = ring::fromInteger<T>(x);
        return f(x_prime) >= 0 && (f(x_prime + 1) < 0 || x_prime - 1 < tm);
    };
    Integer t1_prime = ring::floor_of(tm) + rootradix_prime;
    Integer t1;
    if (is_solution1(t1_prime + 2))
    {
        t1 = t1_prime + 2;
    }
    else if (is_solution1(t1_prime + 1))
    {
        t1 = t1_prime + 1;
    }
    else
    {
        t1 = t1_prime;
    }
    Integer t0_prime = ring::ceiling_of(tm) - rootradix_prime;
    Integer t0;
    if (is_solution0(t0_prime - 2))
    {
        t0 = t0_prime - 2;
    }
    else if (is_solution0(t0_prime - 1))
    {
        t0 = t0_prime - 1;
    }
    else
    {
        t0 = t0_prime;
    }
    return Maybe<Pair<Integer>>(Pair<Integer>{t0, t1});
}

template <typename T, typename U>
Maybe<Pair<U>> quadratic_fixedprec(T a, T b, T c)
{
    int d = std::numeric_limits<U>::digits10;
    std::cout << "d = " << d << std::endl;
    Integer prec = ring::powNonNeg<Integer>(10, d);
    std::cout << "prec = " << prec << std::endl;
    T b_prime = T(prec) * (b * ring::recip(a));
    T c_prime = T(prec) * T(prec) * (c * ring::recip(a));
    std::cout << "b_prime = " << b_prime << std::endl;
    std::cout << "c_prime = " << c_prime << std::endl;
    Maybe<Pair<Integer>> q = int_quadratic(b_prime, c_prime);
    if (!q.has_value())
    {
        return Maybe<Pair<U>>();
    }
    std::cout << "quadratic_fixedprec: q.value = " << q.value() << std::endl;
    Integer x0, x1;
    std::tie(x0, x1) = q.value();
    U prec_u = ring::fromInteger<U>(prec);
    U res1 = ring::fromInteger<U>(x0) * ring::recip(prec_u);
    U res2 = ring::fromInteger<U>(x1) * ring::recip(prec_u);
    return Maybe<Pair<U>>(Pair<U>{res1, res2});
}

template <typename T, typename U>
Maybe<Pair<U>> quadratic(T a, T b, T c)
{
    return quadratic_fixedprec<T, U>(a, b, c);
}