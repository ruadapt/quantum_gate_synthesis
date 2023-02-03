#include "euclideanDomain.h"
#include "ring.h"
#include <tuple>

namespace euclidean_domain
{
    Integer rank(Integer x)
    {
        return x;
    }

    Integer rank(ZOmega x)
    {
        return abs(ring::norm(x));
    }

    Integer rank(ZComplex x)
    {
        return abs(ring::norm(x));
    }

    Integer rank(ZRootTwo x)
    {
        return abs(ring::norm(x));
    }

    std::tuple<Integer, Integer> divmod(Integer x, Integer y)
    {
        return divMod(x, y);
    }

    std::tuple<ZOmega, ZOmega> divmod(ZOmega x, ZOmega y)
    {
        ZOmega o = x * y.adj() * (y * y.adj()).adj2();
        Integer k = ring::norm(y);
        Integer a = rounddiv(o.a(), k);
        Integer b = rounddiv(o.b(), k);
        Integer c = rounddiv(o.c(), k);
        Integer d = rounddiv(o.d(), k);
        ZOmega q = ZOmega(a, b, c, d);
        ZOmega r = x - y * q;
        return std::make_tuple(q, r);
    }

    std::tuple<ZComplex, ZComplex> divmod(ZComplex x, ZComplex y)
    {
        ZComplex c = x * y.adj();
        Integer l = c.a();
        Integer m = c.b();
        Integer k = ring::norm(y);
        Integer q1 = rounddiv(l, k);
        Integer q2 = rounddiv(m, k);
        ZComplex q = ZComplex(q1, q2);
        ZComplex r = x - y * q;
        return std::make_tuple(q, r);
    }

    std::tuple<ZRootTwo, ZRootTwo> divmod(ZRootTwo x, ZRootTwo y)
    {
        ZRootTwo temp = x * y.adj2();
        Integer l = temp.a();
        Integer m = temp.b();
        Integer k = ring::norm(y);
        Integer q1 = rounddiv(l, k);
        Integer q2 = rounddiv(m, k);
        ZRootTwo q = ZRootTwo(q1, q2);
        ZRootTwo r = x - y * q;
        return std::make_tuple(q, r);
    }

    template <typename T>
    T euclid_mod(T x, T y)
    {
        return std::get<1>(divmod(x, y));
    }

    template <typename T>
    T euclid_div(T x, T y)
    {
        return std::get<0>(divmod(x, y));
    }

    template <typename T>
    T euclid_gcd(T x, T y)
    {
        if (y == 0)
        {
            return x;
        }
        T r = std::get<1>(divmod(x, y));
        return euclid_gcd(y, r);
    }

    template <typename T>
    std::tuple<T, T, T, T, T> extended_euclid(T x, T y)
    {
        if (y == 0)
        {
            return std::make_tuple(T(1), T(0), T(0), T(1), x);
        }
        T q, r;
        std::tie(q, r) = divmod(x, y);
        T a2, b2, s2, t2, d;
        std::tie(a2, b2, s2, t2, d) = extended_euclid(y, r);
        return std::make_tuple(b2, a2 - b2 * q, -t2, t2 * q - s2, d);
    }

    template <typename T>
    Maybe<T> euclid_inverse(T x)
    {
        if (x == 0)
        {
            return std::nullopt;
        }
        T q, r;
        std::tie(q, r) = divmod(T(1), x);
        if (r == 0)
        {
            return q;
        }
        return std::nullopt;
    }

    template <typename T>
    bool is_unit(T x)
    {
        return euclid_inverse(x).has_value();
    }

    template <typename T>
    Maybe<T> inv_mod(T p, T a)
    {
        T b, d;
        std::tie(b, std::ignore, std::ignore, std::ignore, d) = extended_euclid(a, p);
        Maybe<T> inv = euclid_inverse(d);
        if (inv.has_value())
        {
            T d2 = inv.value();
            return std::get<1>(divmod(b * d2, p));
        }
        return std::nullopt;
    }

    template <typename T>
    bool euclid_divides(T a, T b)
    {
        if (a == T(0) && b == T(0))
        {
            return true;
        }
        if (a == T(0))
        {
            return false;
        }
        return euclid_mod(b, a) == T(0);
    }

    template <typename T>
    bool euclid_associates(T a, T b)
    {
        return euclid_divides(a, b) && euclid_divides(b, a);
    }

    template <typename T>
    std::tuple<Integer, T> euclid_extract_power(T x, T y)
    {
        if (x == T(0))
        {
            return std::make_tuple(0, x);
        }
        if (is_unit(y))
        {
            return std::make_tuple(0, x);
        }
        if (euclid_divides(y, x))
        {
            Integer k;
            T z;
            std::tie(k, z) = euclid_extract_power(euclid_div(x, y), y);
            return std::make_tuple(k + 1, z);
        }
        return std::make_tuple(0, x);
    }

    Integer rounddiv(Integer x, Integer y)
    {
        Integer temp1, temp2, result;
        mpz_tdiv_q(temp1.get_mpz_t(), y.get_mpz_t(), Integer(2).get_mpz_t());
        temp2 = x + temp1;
        mpz_fdiv_q(result.get_mpz_t(), temp2.get_mpz_t(), y.get_mpz_t());
        return result;
    }

    std::tuple<Integer, Integer> divMod(Integer x, Integer y)
    {
        Integer q, r;
        mpz_fdiv_qr(q.get_mpz_t(), r.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
        return std::make_tuple(q, r);
    }
}