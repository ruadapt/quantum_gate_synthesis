#include "ring.h"
#include <iostream>
#include <tuple>
#include <gmpxx.h>

namespace ring
{
    bool even(int n) { return n % 2 == 0; }

    bool even(Integer n) { return n % 2 == 0; }

    bool odd(int n) { return n % 2 != 0; }

    bool odd(Integer n) { return n % 2 != 0; }

    int div(int a, int b)
    {
        if (b == 0)
        {
            throw std::invalid_argument("Division by 0");
        }
        return mpzToInt(div(Integer(a), Integer(b)));
    }

    Integer div(Integer a, Integer b)
    {
        if (b == 0)
        {
            throw std::invalid_argument("Division by 0");
        }
        Integer result;
        // mpz_fdiv_q always rounds down.
        mpz_fdiv_q(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
        return result;
    }

    int mpzToInt(Integer z)
    {
        assert(z.fits_sint_p());
        return static_cast<int>(z.get_si());
    }

    template <typename Integral>
    Integral shift(Integral x, int bits);

    template <>
    int shift(int x, int bits)
    {
        if (bits >= 0)
        {
            return x << bits;
        }
        return x >> -bits;
    }

    template <typename Integral>
    Integral shift(Integral x, Integer bits)
    {
        return shift<Integral>(x, mpzToInt(bits));
    }

    template <>
    Integer shift(Integer x, int bits)
    {
        Integer result;
        if (bits >= 0)
        {
            mpz_mul_2exp(result.get_mpz_t(), x.get_mpz_t(), bits);
        }
        else
        {
            mpz_fdiv_q_2exp(result.get_mpz_t(), x.get_mpz_t(), -bits);
        }
        return result;
    }

    template <>
    Integer shift(Integer x, Integer bits)
    {
        return shift<Integer>(x, mpzToInt(bits));
    }

    template <typename Integral>
    Integral shiftL(Integral x, int bits)
    {
        return shift<Integral>(x, bits);
    }

    template <typename Integral>
    Integral shiftL(Integral x, Integer bits)
    {
        return shift<Integral>(x, bits);
    }

    template <typename Integral>
    Integral shiftR(Integral x, int bits)
    {
        return shift<Integral>(x, -bits);
    }

    template <typename Integral>
    Integral shiftR(Integral x, Integer bits)
    {
        return shift<Integral>(x, -bits);
    }

    template <typename Integral>
    Integral exp2(int pow)
    {
        return shiftL<Integral>(1, pow);
    }

    template <typename Integral>
    Integral exp2(Integer pow)
    {
        return exp2<Integral>(mpzToInt(pow));
    }

    template <typename T>
    int sign(T a)
    {
        return a.signum();
    }

    template <>
    int sign(int a)
    {
        if (a == 0)
        {
            return 0;
        }
        else if (a > 0)
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }

    template <>
    int sign(Integer a)
    {
        if (a == 0)
        {
            return 0;
        }
        else if (a > 0)
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }

    template <>
    int sign(Real a)
    {
        if (a == 0)
        {
            return 0;
        }
        else if (a > 0)
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }

    template <>
    int sign(Rational a)
    {
        return sgn(a);
    }

    bool getBit(int n, int bit)
    {
        return (n >> bit) & 1;
    }

    int lobit(int n)
    {
        if (n == 0)
        {
            return -1;
        }
        int bit = 0;
        while (1)
        {
            if (getBit(n, bit))
            {
                return bit;
            }
            bit++;
        }
    }

    int lobit(Integer n)
    {
        if (n == 0)
        {
            return -1;
        }
        return static_cast<int>(mpz_scan1(n.get_mpz_t(), 0));
    }

    int hibit(int n)
    {
        if (n == 0)
        {
            return -1;
        }
        int bit = (sizeof n) * 8 - 1;
        while (1)
        {
            if (getBit(n, bit))
            {
                return bit;
            }
            bit--;
        }
    }

    size_t hibit(Integer n)
    {
        // We need the - 1 because mpz_sizeinbase returns the position of the first 1
        // bit counting from 1 instead of 0.
        return mpz_sizeinbase(n.get_mpz_t(), 2) - 1;
    }

    /**
     * Return the floor of sqrt(n).
     */
    Integer intsqrt(Integer n)
    {
        if (n <= 0)
        {
            return 0;
        }
        Integer sqrt;
        // This computes the square root truncated to an integer, and since the square
        // root is non-negative, truncating it is equivalent to taking the floor.
        mpz_sqrt(sqrt.get_mpz_t(), n.get_mpz_t());
        return sqrt;
    }

    std::optional<Integer> log2(Integer n)
    {
        if (n <= 0)
        {
            return std::nullopt;
        }
        int k = lobit(n);
        Integer powK = exp2<Integer>(k);
        if (n == powK)
        {
            return k;
        }
        return std::nullopt;
    }

    template <typename T>
    T fromInteger(int arg)
    {
        return T::fromInteger(arg);
    }

    template <>
    int fromInteger(int arg) { return arg; }

    template <>
    Integer fromInteger(int arg) { return Integer(arg); }

    template <>
    double fromInteger(int arg) { return double(arg); }

    template <>
    Real fromInteger(int arg) { return Real(arg); }

    template <>
    Rational fromInteger(int arg) { return Rational(arg); }

    template <typename T>
    T fromInteger(Integer arg)
    {
        return T::fromInteger(arg);
    }

    template <>
    Integer fromInteger(Integer arg) { return arg; }

    template <>
    double fromInteger(Integer arg) { return arg.get_d(); }

    template <>
    Real fromInteger(Integer arg) { return Real(arg.get_d()); }

    template <>
    Rational fromInteger(Integer arg) { return Rational(arg); }

    template <typename T>
    T recip(T arg)
    {
        return arg.recip();
    }

    template <>
    Real recip(Real arg)
    {
        assert(arg != 0);
        return 1 / arg;
    }

    template <>
    Rational recip(Rational arg)
    {
        assert(arg != 0);
        return 1 / arg;
    }

    template <typename T>
    T fromRational(Rational r)
    {
        return T::fromRational(r);
    }

    template <>
    Rational fromRational(Rational r)
    {
        return r;
    }

    template <>
    Real fromRational(Rational r)
    {
        return r.get_d();
    }

    template <typename T>
    T powNonNeg(T base, int exp)
    {
        assert(exp >= 0);
        if (exp == 0)
        {
            return 1;
        }
        if (exp == 1)
        {
            return base;
        }
        if (exp % 2 == 0)
        {
            T newBase = base * base;
            int newExp = exp / 2;
            return powNonNeg(newBase, newExp);
        }
        T newBase = base * base;
        int newExp = (exp - 1) / 2;
        return base * powNonNeg(newBase, newExp);
    }

    template <typename T>
    T powNonNeg(T base, Integer exp)
    {
        // Any exponent should be able to fit in an int, or else the result would
        // be too large.
        return powNonNeg<T>(base, mpzToInt(exp));
    }

    /**
     * This will only work for types that have recip defined.
     */
    template <typename T>
    T powInt(T base, int exp)
    {
        return (exp >= 0) ? powNonNeg<T>(base, exp) : powNonNeg<T>(recip(base), -exp);
    }

    template <typename T>
    T powInt(T base, Integer exp)
    {
        // Any exponent should be able to fit in an int, or else the result would
        // be too large.
        return powInt<T>(base, mpzToInt(exp));
    }

    template <typename T>
    T half()
    {
        return T::half();
    }

    template <>
    Real half() { return 0.5; }

    template <>
    Rational half() { return Rational(1, 2); }

    template <typename T>
    T fromDyadic(Dyadic<int> d)
    {
        std::tuple<int, int> decomposed = d.decomposeDyadic();
        int a = std::get<0>(decomposed);
        int n = std::get<1>(decomposed);
        if (n >= 0)
        {
            return fromInteger<T>(a) * powNonNeg<T>(half<T>(), n);
        }
        else
        {
            return fromInteger<T>(a) * fromInteger<T>(exp2<Integer>(-n));
        }
    }

    template <typename T>
    T fromDyadic(ZDyadic d)
    {
        std::tuple<Integer, Integer> decomposed = d.decomposeDyadic();
        Integer a = std::get<0>(decomposed);
        Integer n = std::get<1>(decomposed);
        if (n >= 0)
        {
            return fromInteger<T>(a) * powNonNeg<T>(half<T>(), n);
        }
        else
        {
            return fromInteger<T>(a) * fromInteger<T>(exp2<Integer>(-n));
        }
    }

    template <typename T>
    T rootTwo()
    {
        return T::rootTwo();
    }

    template <>
    Real rootTwo() { return bmp::sqrt(Real(2)); }

    template <typename T>
    T fromZRootTwo(ZRootTwo arg)
    {
        return fromInteger<T>(arg.a()) + rootTwo<T>() * fromInteger<T>(arg.b());
    }

    template <typename T>
    T rootHalf()
    {
        return T::rootHalf();
    }

    template <>
    Real rootHalf()
    {
        return bmp::sqrt(Real(0.5));
    }

    template <typename T>
    T fromDRootTwo(DRootTwo arg)
    {
        return fromDyadic<T>(arg.a()) + rootTwo<T>() * fromDyadic<T>(arg.b());
    }

    template <typename T>
    T i()
    {
        return T::i();
    }

    template <typename T>
    T adj(T arg)
    {
        return arg.adj();
    }

    template <>
    int adj(int arg) { return arg; }

    template <>
    Integer adj(Integer arg) { return arg; }

    template <>
    Real adj(Real arg) { return arg; }

    template <>
    Rational adj(Rational arg) { return arg; }

    template <typename T>
    T adj2(T arg)
    {
        return arg.adj2();
    }

    template <>
    int adj2(int arg) { return arg; }

    template <>
    Integer adj2(Integer arg) { return arg; }

    template <>
    Real adj2(Real arg) { return arg; }

    template <>
    Rational adj2(Rational arg) { return arg; }

    template <typename T>
    Integer norm(T arg)
    {
        return arg.norm();
    }

    template <>
    Integer norm(int arg) { return mpz_class(arg); }

    template <>
    Integer norm(Integer arg) { return arg; }

    Integer floor_of(double arg)
    {
        return floor(arg);
    }

    Integer floor_of(Real arg)
    {
        return floor(double(arg));
    }

    Integer floor_of(Rational arg)
    {
        return div(arg.get_num(), arg.get_den());
    }

    Integer floor_of(Integer arg)
    {
        return arg;
    }

    Integer floor_of(QRootTwo arg)
    {
        Integer a = floor_of(arg.a());
        Integer b = intsqrt(floor_of(2 * arg.b() * arg.b()));
        Integer rInt = (arg.b() >= 0) ? (Integer(a + b)) : (Integer(a - b));
        QRootTwo r = QRootTwo(Rational(rInt), 0);
        if (r + ring::fromInteger<QRootTwo>(1) <= arg)
        {
            return rInt + 1;
        }
        else if (r <= arg)
        {
            return rInt;
        }
        return rInt - 1;
    }

    Integer ceiling_of(double arg)
    {
        return ceil(arg);
    }

    Integer ceiling_of(Real arg)
    {
        return ceil(double(arg));
    }

    Integer ceiling_of(Rational arg)
    {
        Integer result;
        mpz_cdiv_q(result.get_mpz_t(), arg.get_num_mpz_t(), arg.get_den_mpz_t());
        return result;
    }

    Integer ceiling_of(Integer arg)
    {
        return arg;
    }

    Integer ceiling_of(QRootTwo arg)
    {
        return -floor_of(-arg);
    }

    template <typename T>
    T omega()
    {
        return T::omega();
    }

    template <typename T, typename U>
    std::optional<RootTwo<U>> maybeDyadic(RootTwo<T> arg)
    {
        std::optional<U> aDyadic = ring::maybeDyadic<T, U>(arg.a());
        std::optional<U> bDyadic = ring::maybeDyadic<T, U>(arg.b());
        if (aDyadic.has_value() && bDyadic.has_value())
        {
            return RootTwo<U>(aDyadic.value(), bDyadic.value());
        }
        return std::nullopt;
    }

    template <typename T, typename U>
    std::optional<Complex<U>> maybeDyadic(Complex<T> arg)
    {
        std::optional<U> aDyadic = ring::maybeDyadic<T, U>(arg.a());
        std::optional<U> bDyadic = ring::maybeDyadic<T, U>(arg.b());
        if (aDyadic.has_value() && bDyadic.has_value())
        {
            return Complex<U>(aDyadic.value(), bDyadic.value());
        }
        return std::nullopt;
    }

    template <typename T, typename U>
    std::optional<Omega<U>> maybeDyadic(Omega<T> arg)
    {
        std::optional<U> aDyadic = ring::maybeDyadic<T, U>(arg.a());
        std::optional<U> bDyadic = ring::maybeDyadic<T, U>(arg.b());
        std::optional<U> cDyadic = ring::maybeDyadic<T, U>(arg.c());
        std::optional<U> dDyadic = ring::maybeDyadic<T, U>(arg.d());
        if (aDyadic.has_value() && bDyadic.has_value() && cDyadic.has_value() && dDyadic.has_value())
        {
            return Omega<U>(aDyadic.value(), bDyadic.value(), cDyadic.value(), dDyadic.value());
        }
        return std::nullopt;
    }

    template <>
    std::optional<ZDyadic> maybeDyadic(ZDyadic arg)
    {
        return arg;
    }

    template <>
    std::optional<ZDyadic> maybeDyadic(Rational r)
    {
        std::optional<Integer> k = log2(r.get_den());
        if (k.has_value())
        {
            return ZDyadic(r.get_num(), k.value());
        }
        return std::nullopt;
    }

    template <typename T, typename U>
    U toDyadic(T arg)
    {
        std::optional<U> maybe = maybeDyadic<T, U>(arg);
        if (maybe.has_value())
        {
            return maybe.value();
        }
        throw std::invalid_argument("Could not convert argument to dyadic value");
    }

    template <typename T, typename U>
    RootTwo<U> toDyadic(RootTwo<T> arg)
    {
        std::optional<RootTwo<U>> maybe = maybeDyadic<T, U>(arg);
        if (maybe.has_value())
        {
            return maybe.value();
        }
        throw std::invalid_argument("Could not convert argument to dyadic value");
    }

    template <typename T, typename U>
    Complex<U> toDyadic(Complex<T> arg)
    {
        std::optional<Complex<U>> maybe = maybeDyadic<T, U>(arg);
        if (maybe.has_value())
        {
            return maybe.value();
        }
        throw std::invalid_argument("Could not convert argument to dyadic value");
    }

    template <typename T, typename U>
    Omega<U> toDyadic(Omega<T> arg)
    {
        std::optional<Omega<U>> maybe = maybeDyadic<T, U>(arg);
        if (maybe.has_value())
        {
            return maybe.value();
        }
        throw std::invalid_argument("Could not convert argument to dyadic value");
    }

    template <typename T>
    T real(Complex<T> arg)
    {
        return arg.a();
    }

    template <typename T>
    RootTwo<T> real(Omega<T> arg)
    {
        return RootTwo<T>(arg.d(), half<T>() * (arg.c() - arg.a()));
    }

    template <typename T, typename U>
    T fromWhole(U arg);

    template <>
    ZDyadic fromWhole(Integer arg)
    {
        return fromInteger<ZDyadic>(arg);
    }

    template <>
    DOmega fromWhole(ZOmega arg)
    {
        return fromZOmega<DOmega>(arg);
    }

    template <>
    DRootTwo fromWhole(ZRootTwo arg)
    {
        return fromZRootTwo<DRootTwo>(arg);
    }

    template <typename T, typename U>
    U toWhole(T arg);

    template <>
    Integer toWhole(ZDyadic arg)
    {
        Integer a, n;
        std::tie(a, n) = arg.decomposeDyadic();
        if (n == 0)
        {
            return a;
        }
        throw std::invalid_argument("Non-integral value can't be converted to Integer");
    }

    template <>
    ZOmega toWhole(DOmega arg)
    {
        Integer a2 = toWhole<ZDyadic, Integer>(arg.a());
        Integer b2 = toWhole<ZDyadic, Integer>(arg.b());
        Integer c2 = toWhole<ZDyadic, Integer>(arg.c());
        Integer d2 = toWhole<ZDyadic, Integer>(arg.d());
        return ZOmega(a2, b2, c2, d2);
    }

    template <>
    ZRootTwo toWhole(DRootTwo arg)
    {
        Integer a2 = toWhole<ZDyadic, Integer>(arg.a());
        Integer b2 = toWhole<ZDyadic, Integer>(arg.b());
        return ZRootTwo(a2, b2);
    }

    template <typename T>
    Integer denomExp(T arg)
    {
        return arg.denomExp();
    }

    template <>
    Integer denomExp(DOmega arg)
    {
        Integer a, ak, b, bk, c, ck, d, dk;
        std::tie(a, ak) = arg.a().decomposeDyadic();
        std::tie(b, bk) = arg.b().decomposeDyadic();
        std::tie(c, ck) = arg.c().decomposeDyadic();
        std::tie(d, dk) = arg.d().decomposeDyadic();
        Integer maxK = std::max({ak, bk, ck, ck});
        Integer a2 = (maxK == ak) ? a : 0;
        Integer b2 = (maxK == bk) ? b : 0;
        Integer c2 = (maxK == ck) ? c : 0;
        Integer d2 = (maxK == dk) ? d : 0;
        Integer k2;
        if (maxK > 0 && even(a2 - c2) && even(b2 - d2))
        {
            return 2 * maxK - 1;
        }
        return 2 * maxK;
    }

    template <>
    Integer denomExp(DRootTwo arg)
    {
        Integer a, ak, b, bk;
        std::tie(a, ak) = arg.a().decomposeDyadic();
        std::tie(b, bk) = arg.b().decomposeDyadic();
        Integer x1 = 2 * ak;
        Integer x2 = 2 * bk - 1;
        return std::max<Integer>(x1, x2);
    }

    template <typename T>
    T denomExpFactor(T arg, Integer k)
    {
        return arg.denomExpFactor(k);
    }

    template <>
    DOmega denomExpFactor(DOmega arg, Integer k)
    {
        return arg * powNonNeg(rootTwo<DOmega>(), k);
    }

    template <>
    DRootTwo denomExpFactor(DRootTwo arg, Integer k)
    {
        return arg * powNonNeg(rootTwo<DRootTwo>(), k);
    }

    template <typename T>
    QOmega toQOmega(T arg)
    {
        return arg.toQOmega();
    }

    template <>
    QOmega toQOmega(int arg) { return fromInteger<QOmega>(arg); }

    template <>
    QOmega toQOmega(Integer arg) { return fromInteger<QOmega>(arg); }

    template <>
    QOmega toQOmega(Rational arg) { return fromRational<QOmega>(arg); }

    Z2 parity(int arg)
    {
        return (arg % 2) != 0;
    }

    Z2 parity(Integer arg)
    {
        return (arg % 2) != 0;
    }

    Z2 parity(ZRootTwo arg)
    {
        return parity(arg.a());
    }

    template <typename T>
    T fromQRootTwo(QRootTwo q)
    {
        return fromRational<T>(q.a()) + rootTwo<T>() * fromRational<T>(q.b());
    }

    template <typename T>
    T fromZComplex(ZComplex z)
    {
        return fromInteger<T>(z.a()) + i<T>() * fromInteger<T>(z.b());
    }

    template <typename T>
    T fromDComplex(DComplex d)
    {
        return fromDyadic<T>(d.a()) + i<T>() * fromDyadic<T>(d.b());
    }

    template <typename T>
    T fromQComplex(QComplex q)
    {
        return fromRational<T>(q.a()) + i<T>() * fromRational<T>(q.b());
    }

    template <typename T>
    T fromDRComplex(DRComplex d)
    {
        return fromDRootTwo<T>(d.a()) + i<T>() * fromDRootTwo<T>(d.b());
    }

    template <typename T>
    T fromQRComplex(QRComplex q)
    {
        return fromQRootTwo<T>(q.a()) + i<T>() * fromQRootTwo<T>(q.b());
    }

    template <typename T>
    T fromZOmega(ZOmega z)
    {
        T o = omega<T>();
        T o2 = o * o;
        T o3 = o2 * o;
        return fromInteger<T>(z.a()) * o3 + fromInteger<T>(z.b()) * o2 + fromInteger<T>(z.c()) * o + fromInteger<T>(z.d());
    }

    template <typename T>
    T fromDOmega(DOmega z)
    {
        T o = omega<T>();
        T o2 = o * o;
        T o3 = o2 * o;
        return fromDyadic<T>(z.a()) * o3 + fromDyadic<T>(z.b()) * o2 + fromDyadic<T>(z.c()) * o + fromDyadic<T>(z.d());
    }

    template <typename T>
    T fromQOmega(QOmega z)
    {
        T o = omega<T>();
        T o2 = o * o;
        T o3 = o2 * o;
        return fromRational<T>(z.a()) * o3 + fromRational<T>(z.b()) * o2 + fromRational<T>(z.c()) * o + fromRational<T>(z.d());
    }

    std::optional<ZRootTwo> zRootTwoRoot(ZRootTwo arg)
    {
        Integer d = arg.a() * arg.a() - 2 * arg.b() * arg.b();
        Integer r = intsqrt(d);
        Integer sum = arg.a() + r;
        Integer diff = arg.a() - r;
        Integer integer2 = 2;
        Integer integer4 = 4;
        Integer div1 = div(sum, integer2);
        Integer div2 = div(diff, integer2);
        Integer div3 = div(diff, integer4);
        Integer div4 = div(sum, integer4);
        Integer x1 = intsqrt(div1);
        Integer x2 = intsqrt(div2);
        Integer y1 = intsqrt(div3);
        Integer y2 = intsqrt(div4);
        ZRootTwo w1 = ZRootTwo(x1, y1);
        ZRootTwo w2 = ZRootTwo(x2, y2);
        ZRootTwo w3 = ZRootTwo(x1, -y1);
        ZRootTwo w4 = ZRootTwo(x2, -y2);
        if (w1 * w1 == arg)
        {
            return w1;
        }
        if (w2 * w2 == arg)
        {
            return w2;
        }
        if (w3 * w3 == arg)
        {
            return w3;
        }
        if (w4 * w4 == arg)
        {
            return w4;
        }
        return std::nullopt;
    }

    ZRootTwo zRootTwoOfZOmega(ZOmega arg)
    {
        if ((arg.a() == -arg.c()) && (arg.b() == 0))
        {
            return ZRootTwo(arg.d(), arg.c());
        }
        throw std::invalid_argument("Could not convert ZOmega value to ZRootTwo");
    }

    template <typename T>
    std::string toString(const T &arg)
    {
        return arg.toString();
    }

    template <>
    std::string toString(const int &arg) { return std::to_string(arg); }

    template <>
    std::string toString(const Integer &arg) { return arg.get_str(); }

    template <>
    std::string toString(const double &arg) { return std::to_string(arg); }

    std::string toString(const Real &arg)
    {
        std::stringstream ss;
        ss << std::setprecision(std::numeric_limits<Real>::digits10) << arg;
        return ss.str();
    }

    template <>
    std::string toString(const Rational &arg) { return arg.get_str(); }
}

template <typename T>
Dyadic<T>::Dyadic()
{
    a_ = 0;
    n_ = 0;
}

template <typename T>
Dyadic<T>::Dyadic(int arg)
{
    a_ = arg;
    n_ = 0;
}

template <>
ZDyadic::Dyadic(Integer arg)
{
    a_ = arg;
    n_ = 0;
}

template <typename T>
Dyadic<T> &Dyadic<T>::operator+=(const Dyadic<T> &d)
{
    Dyadic<T> sum = *this + d;
    this->a = sum.a;
    this->n = sum.n;
    return *this;
}

template <typename T>
Dyadic<T>::Dyadic(T a, T n)
{
    a_ = a;
    n_ = n;
}

template <typename T>
T Dyadic<T>::a() const
{
    return a_;
}

template <typename T>
T Dyadic<T>::n() const
{
    return n_;
}

template <typename T>
Dyadic<T> Dyadic<T>::copy() const
{
    return Dyadic<T>(a(), n());
}

template <typename T>
bool Dyadic<T>::operator==(const Dyadic &d) const
{
    T b = d.a();
    T m = d.n();
    T k = std::max<T>(m, n());
    return (a() * ring::exp2<T>(k - n())) == (b * ring::exp2<T>(k - m));
}

template <typename T>
bool Dyadic<T>::operator!=(const Dyadic &d) const
{
    return !(*this == d);
}

template <typename T>
Ordering Dyadic<T>::compare(const Dyadic &d) const
{
    T b = d.a();
    T m = d.n();
    T k = std::max<T>(n(), m);
    T size1 = a() * ring::exp2<T>(k - n());
    T size2 = b * ring::exp2<T>(k - m);
    if (size1 > size2)
    {
        return GT;
    }
    else if (size1 < size2)
    {
        return LT;
    }
    return EQ;
}

template <typename T>
bool Dyadic<T>::operator<(const Dyadic &d) const
{
    return compare(d) == LT;
}

template <typename T>
bool Dyadic<T>::operator>(const Dyadic &d) const
{
    return compare(d) == GT;
}

template <typename T>
bool Dyadic<T>::operator<=(const Dyadic &d) const
{
    Ordering c = compare(d);
    return (c == LT) || (c == EQ);
}

template <typename T>
bool Dyadic<T>::operator>=(const Dyadic &d) const
{
    Ordering c = compare(d);
    return (c == GT) || (c == EQ);
}

template <typename T>
Dyadic<T> Dyadic<T>::operator+(const Dyadic &d) const
{
    T b = d.a();
    T m = d.n();
    if (n() < m)
    {
        T c = ring::shiftL(a(), m - n()) + b;
        return Dyadic(c, m);
    }
    else
    {
        T d = a() + ring::shiftL(b, n() - m);
        return Dyadic(d, n());
    }
}

template <typename T>
Dyadic<T> Dyadic<T>::operator-(const Dyadic &d) const
{
    return (*this) + (-d);
}

template <typename T>
Dyadic<T> Dyadic<T>::operator*(const Dyadic &d) const
{
    T b = d.a();
    T m = d.n();
    return Dyadic(a() * b, m + n());
}

template <typename T>
Dyadic<T> Dyadic<T>::operator-() const
{
    return Dyadic(-a(), n());
}

template <typename T>
Dyadic<T> Dyadic<T>::abs() const
{
    if (compare(ring::fromInteger<Dyadic<T>>(0)) != LT)
    {
        return copy();
    }
    return -copy();
}

template <typename T>
int Dyadic<T>::signum() const
{
    Ordering comp = compare(ring::fromInteger<Dyadic<T>>(0));
    if (comp == LT)
    {
        return -1;
    }
    else if (comp == GT)
    {
        return 1;
    }
    return 0;
}

template <typename T>
Dyadic<T> Dyadic<T>::adj() const
{
    return copy();
}

template <typename T>
Dyadic<T> Dyadic<T>::adj2() const
{
    return copy();
}

template <typename T>
std::tuple<T, T> Dyadic<T>::decomposeDyadic() const
{
    if (a() == 0)
    {
        return std::make_tuple(0, 0);
    }
    else
    {
        int k = ring::lobit(a());
        if (n() >= k)
        {
            return std::make_tuple(ring::shiftR(a(), k), n() - k);
        }
        else
        {
            return std::make_tuple(ring::shiftR(a(), n()), 0);
        }
    }
}

template <typename T>
T Dyadic<T>::integerOfDyadic(T k) const
{
    return ring::shift(a(), (k - n()));
}

template <typename T>
QOmega Dyadic<T>::toQOmega() const
{
    if (n() >= 0)
    {
        return ring::toQOmega<T>(a()) * ring::powNonNeg<QOmega>(ring::half<QOmega>(), n());
    }
    return ring::toQOmega<T>(a()) * ring::powNonNeg<QOmega>(QOmega(2), -n());
}

template <typename T>
std::string Dyadic<T>::toString() const
{
    return "Dyadic(" + ring::toString(a()) + ", " + ring::toString(n()) + ")";
}

template <typename T>
void Dyadic<T>::print(std::string prefix) const
{
    std::cout << prefix << ": " << toString() << std::endl;
}

template <typename T>
Dyadic<T> Dyadic<T>::fromInteger(int n)
{
    return Dyadic(n, 0);
}

template <typename T>
Dyadic<T> Dyadic<T>::fromInteger(Integer n)
{
    return Dyadic(ring::fromInteger<T>(n), 0);
}

template <typename T>
Dyadic<T> Dyadic<T>::fromDyadic(const Dyadic &d)
{
    return d.copy();
}

template <typename T>
Dyadic<T> Dyadic<T>::half()
{
    return Dyadic(1, 1);
}

template <>
Dyadic<int> ring::half()
{
    return Dyadic<int>(1, 1);
}

template <>
ZDyadic ring::half()
{
    return ZDyadic(1, 1);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Dyadic<T> &d)
{
    os << d.toString();
    return os;
}

template <typename T>
RootTwo<T>::RootTwo()
{
    a_ = 0;
    b_ = 0;
}

template <typename T>
RootTwo<T>::RootTwo(int arg)
{
    a_ = arg;
    b_ = 0;
}

template <typename T>
RootTwo<T>::RootTwo(Integer arg)
{
    a_ = arg;
    b_ = 0;
}

template <typename T>
RootTwo<T>::RootTwo(T a, T b)
{
    a_ = a;
    b_ = b;
}

template <typename T>
RootTwo<T> &RootTwo<T>::operator+=(const RootTwo<T> &r)
{
    RootTwo<T> sum = *this + r;
    this->a_ = sum.a();
    this->b_ = sum.b();
    return *this;
}

template <typename T>
T RootTwo<T>::a() const
{
    return a_;
}

template <typename T>
T RootTwo<T>::b() const
{
    return b_;
}

template <typename T>
RootTwo<T> RootTwo<T>::copy() const
{
    return RootTwo(a(), b());
}

template <typename T>
bool RootTwo<T>::operator==(const RootTwo &r) const
{
    return (a() == r.a()) && (b() == r.b());
}

template <typename T>
bool RootTwo<T>::operator!=(const RootTwo &r) const
{
    return !(*this == r);
}

template <typename T>
bool RootTwo<T>::operator<=(const RootTwo &r) const
{
    return (r - *this).signum() != -1;
}

template <typename T>
bool RootTwo<T>::operator<(const RootTwo &r) const
{
    return (*this <= r) && !(*this == r);
}

template <typename T>
bool RootTwo<T>::operator>=(const RootTwo &r) const
{
    return !(*this < r);
}

template <typename T>
bool RootTwo<T>::operator>(const RootTwo &r) const
{
    return !(*this <= r);
}

template <typename T>
RootTwo<T> RootTwo<T>::operator+(const RootTwo &r) const
{
    return RootTwo(a() + r.a(), b() + r.b());
}

template <typename T>
RootTwo<T> RootTwo<T>::operator-(const RootTwo &r) const
{
    return RootTwo(a() - r.a(), b() - r.b());
}

template <typename T>
RootTwo<T> RootTwo<T>::operator*(const RootTwo &r) const
{
    T newA = a() * r.a() + (b() * r.b()) + (b() * r.b());
    T newB = a() * r.b() + r.a() * b();
    return RootTwo(newA, newB);
}

template <typename T>
RootTwo<T> RootTwo<T>::operator-() const
{
    return RootTwo(-a(), -b());
}

template <typename T>
RootTwo<T> RootTwo<T>::abs() const
{
    int sign = signum();
    if (sign != -1)
    {
        return copy();
    }
    return -(*this);
}

template <typename T>
int RootTwo<T>::signum() const
{
    int sa = ring::sign(a());
    int sb = ring::sign(b());
    if (sa == 0 && sb == 0)
    {
        return 0;
    }
    else if (sa != -1 && sb != -1)
    {
        return 1;
    }
    else if (sa != 1 && sb != 1)
    {
        return -1;
    }
    else if (sa != -1 && sb != 1 && ring::sign<T>(a() * a() - T(2) * b() * b()) != -1)
    {
        return 1;
    }
    else if (sa != 1 && sb != -1 && ring::sign<T>(a() * a() - T(2) * b() * b()) != 1)
    {
        return 1;
    }
    return -1;
}

template <typename T>
RootTwo<T> RootTwo<T>::adj() const
{
    return RootTwo(ring::adj(a()), ring::adj(b()));
}

template <typename T>
RootTwo<T> RootTwo<T>::adj2() const
{
    return RootTwo(ring::adj2(a()), -ring::adj2(b()));
}

template <typename T>
RootTwo<T> RootTwo<T>::recip() const
{
    assert((*this) != 0);
    T k = a() * a() - T(2) * b() * b();
    return RootTwo(a() * ring::recip(k), -b() * ring::recip(k));
}

template <typename T>
Integer RootTwo<T>::norm() const
{
    Integer normA = ring::norm<T>(a());
    Integer normB = ring::norm<T>(b());
    return normA * normA - 2 * normB * normB;
}

template <typename T>
QOmega RootTwo<T>::toQOmega() const
{
    return ring::toQOmega<T>(a()) + ring::rootTwo<QOmega>() * ring::toQOmega<T>(b());
}

template <typename T>
std::string RootTwo<T>::toString() const
{
    return "RootTwo(" + ring::toString(a()) + ", " + ring::toString(b()) + ")";
}

template <typename T>
void RootTwo<T>::print(std::string prefix) const
{
    std::cout << prefix << ": " << toString() << std::endl;
}

template <typename T>
RootTwo<T> RootTwo<T>::half()
{
    return RootTwo<T>(ring::half<T>(), 0);
}

template <typename T>
RootTwo<T> RootTwo<T>::rootTwo()
{
    return RootTwo<T>(0, 1);
}

template <typename T>
RootTwo<T> RootTwo<T>::rootHalf()
{
    return RootTwo<T>(0, ring::half<T>());
}

template <typename T>
RootTwo<T> RootTwo<T>::i()
{
    return RootTwo<T>(ring::i<T>(), 0);
}

template <typename T>
RootTwo<T> RootTwo<T>::omega()
{
    return RootTwo<T>::rootHalf() * (RootTwo<T>(1) + RootTwo<T>::i());
}

template <typename T>
RootTwo<T> RootTwo<T>::fromInteger(int n)
{
    return RootTwo<T>(n, 0);
}

template <typename T>
RootTwo<T> RootTwo<T>::fromInteger(Integer n)
{
    return RootTwo<T>(ring::fromInteger<T>(n), 0);
}

template <typename T>
RootTwo<T> RootTwo<T>::fromRational(Rational r)
{
    return RootTwo<T>(ring::fromRational<T>(r), 0);
}

template <>
QRootTwo::RootTwo(Rational a, Rational b)
{
    // Make sure Rationals are in canonical form.
    a.canonicalize();
    b.canonicalize();
    a_ = a;
    b_ = b;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const RootTwo<T> &r)
{
    os << r.toString();
    return os;
}

template <typename T>
Complex<T>::Complex()
{
    a_ = 0;
    b_ = 0;
}

template <typename T>
Complex<T>::Complex(int arg)
{
    a_ = arg;
    b_ = 0;
}

template <typename T>
Complex<T>::Complex(Integer arg)
{
    a_ = arg;
    b_ = 0;
}

template <typename T>
Complex<T>::Complex(T a, T b)
{
    a_ = a;
    b_ = b;
}

template <typename T>
Complex<T> &Complex<T>::operator+=(const Complex<T> &c)
{
    Complex<T> sum = *this + c;
    this->a_ = sum.a();
    this->b_ = sum.b();
    return *this;
}

template <typename T>
T Complex<T>::a() const
{
    return a_;
}

template <typename T>
T Complex<T>::b() const
{
    return b_;
}

template <typename T>
Complex<T> Complex<T>::copy() const
{
    return Complex(a(), b());
}

template <typename T>
bool Complex<T>::operator==(const Complex &c) const
{
    return (a() == c.a()) && (b() == c.b());
}

template <typename T>
bool Complex<T>::operator!=(const Complex &c) const
{
    return !(*this == c);
}

template <typename T>
Complex<T> Complex<T>::operator+(const Complex &c) const
{
    return Complex(a() + c.a(), b() + c.b());
}

template <typename T>
Complex<T> Complex<T>::operator-(const Complex &c) const
{
    return Complex(a() - c.a(), b() - c.b());
}

template <typename T>
Complex<T> Complex<T>::operator*(const Complex &c) const
{
    T newA = a() * c.a() - b() * c.b();
    T newB = a() * c.b() + c.a() * b();
    return Complex(newA, newB);
}

template <typename T>
Complex<T> Complex<T>::operator-() const
{
    return Complex(-a(), -b());
}

template <typename T>
Complex<T> Complex<T>::abs() const
{
    return copy();
}

template <typename T>
int Complex<T>::signum() const
{
    return 1;
}

template <typename T>
Complex<T> Complex<T>::adj() const
{
    return Complex(ring::adj(a()), -ring::adj(b()));
}

template <typename T>
Complex<T> Complex<T>::adj2() const
{
    return Complex(ring::adj2(a()), ring::adj2(b()));
}

template <typename T>
Complex<T> Complex<T>::recip() const
{
    assert((*this) != 0);
    T d = a() * a() + b() * b();
    return Complex(a() * ring::recip(d), -b() * ring::recip(d));
}

template <typename T>
Integer Complex<T>::norm() const
{
    Integer normA = ring::norm<T>(a());
    Integer normB = ring::norm<T>(b());
    return normA * normA + normB * normB;
}

template <typename T>
Integer Complex<T>::denomExp() const
{
    Integer expA = ring::denomExp<T>(a());
    Integer expB = ring::denomExp<T>(b());
    return std::max<Integer>(expA, expB);
}

template <typename T>
Complex<T> Complex<T>::denomExpFactor(Integer k) const
{
    return Complex<T>(ring::denomExpFactor<T>(a(), k), ring::denomExpFactor<T>(b(), k));
}

template <typename T>
QOmega Complex<T>::toQOmega() const
{
    return ring::toQOmega(a()) + ring::i<QOmega>() * ring::toQOmega<T>(b());
}

template <typename T>
std::string Complex<T>::toString() const
{
    return "Complex(" + ring::toString(a()) + ", " + ring::toString(b()) + ")";
}

template <typename T>
void Complex<T>::print(std::string prefix) const
{
    std::cout << prefix << ": " << toString() << std::endl;
}

template <typename T>
Complex<T> Complex<T>::half()
{
    return Complex<T>(ring::half<T>(), 0);
}

template <typename T>
Complex<T> Complex<T>::rootTwo()
{
    return Complex<T>(ring::rootTwo<T>(), 0);
}

template <typename T>
Complex<T> Complex<T>::rootHalf()
{
    return Complex<T>(ring::rootHalf<T>(), 0);
}

template <typename T>
Complex<T> Complex<T>::i()
{
    return Complex<T>(0, 1);
}

template <typename T>
Complex<T> Complex<T>::omega()
{
    return Complex<T>::rootHalf() * (Complex<T>(1) + Complex<T>::i());
}

template <typename T>
Complex<T> Complex<T>::fromInteger(int n)
{
    return Complex<T>(n, 0);
}

template <typename T>
Complex<T> Complex<T>::fromInteger(Integer n)
{
    return Complex<T>(ring::fromInteger<T>(n), 0);
}

template <typename T>
Complex<T> Complex<T>::fromRational(Rational r)
{
    return Complex<T>(ring::fromRational<T>(r), 0);
}

template <>
QComplex::Complex(Rational a, Rational b)
{
    // Make sure Rationals are in canonical form.
    a.canonicalize();
    b.canonicalize();
    a_ = a;
    b_ = b;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Complex<T> &c)
{
    os << c.toString();
    return os;
}

Z2::Z2()
{
    mod2_ = 0;
}

Z2::Z2(int arg)
{
    mod2_ = (arg % 2) != 0;
}

Z2::Z2(Integer arg)
{
    mod2_ = mpz_odd_p(arg.get_mpz_t());
}

Z2::Z2(bool mod2)
{
    mod2_ = mod2;
}

bool Z2::mod2() const
{
    return mod2_;
}

Z2 Z2::copy() const
{
    return Z2(mod2());
}

bool Z2::operator==(const Z2 &z) const
{
    return mod2() == z.mod2();
}

bool Z2::operator!=(const Z2 &z) const
{
    return !(*this == z);
}

Z2 Z2::operator+(const Z2 &z) const
{
    return Z2(mod2() != z.mod2());
}

Z2 Z2::operator-(const Z2 &z) const
{
    return (*this) + (-z);
}

Z2 Z2::operator*(const Z2 &z) const
{
    if (mod2() == 0)
    {
        return Z2(0);
    }
    return z.mod2();
}

Z2 Z2::operator-() const
{
    return copy();
}
Z2 Z2::abs() const
{
    return copy();
}

int Z2::signum() const
{
    return 1;
}

Z2 Z2::adj() const
{
    return copy();
}

Z2 Z2::adj2() const
{
    return copy();
}

std::string Z2::toString() const
{
    return "Z2(" + std::to_string(mod2()) + ")";
}

void Z2::print(std::string prefix) const
{
    std::cout << prefix << ": " << toString() << std::endl;
}

Z2 Z2::fromInteger(int n)
{
    return Z2(n);
}

Z2 Z2::fromInteger(Integer n)
{
    return Z2(n);
}

std::ostream& operator<<(std::ostream& os, const Z2 &z)
{
    os << z.toString();
    return os;
}

template <typename T>
Omega<T>::Omega()
{
    a_ = 0;
    b_ = 0;
    c_ = 0;
    d_ = 0;
}

template <typename T>
Omega<T>::Omega(int arg)
{
    a_ = 0;
    b_ = 0;
    c_ = 0;
    d_ = arg;
}

template <typename T>
Omega<T>::Omega(Integer arg)
{
    a_ = 0;
    b_ = 0;
    c_ = 0;
    d_ = arg;
}

template <typename T>
Omega<T>::Omega(T a, T b, T c, T d)
{
    a_ = a;
    b_ = b;
    c_ = c;
    d_ = d;
}

template <typename T>
Omega<T> &Omega<T>::operator+=(const Omega<T> &o)
{
    Omega<T> sum = *this + o;
    this->a_ = sum.a();
    this->b_ = sum.b();
    this->c_ = sum.c();
    this->d_ = sum.d();
    return *this;
}

template <typename T>
T Omega<T>::a() const
{
    return a_;
}

template <typename T>
T Omega<T>::b() const
{
    return b_;
}

template <typename T>
T Omega<T>::c() const
{
    return c_;
}

template <typename T>
T Omega<T>::d() const
{
    return d_;
}

template <typename T>
Omega<T> Omega<T>::copy() const
{
    return Omega<T>(a(), b(), c(), d());
}

template <typename T>
bool Omega<T>::operator==(const Omega &o) const
{
    return (a() == o.a()) && (b() == o.b()) && (c() == o.c()) && (d() == o.d());
}

template <typename T>
bool Omega<T>::operator!=(const Omega &o) const
{
    return !(*this == o);
}

template <typename T>
Omega<T> Omega<T>::operator+(const Omega &o) const
{
    return Omega<T>(a() + o.a(), b() + o.b(), c() + o.c(), d() + o.d());
}

template <typename T>
Omega<T> Omega<T>::operator-(const Omega &o) const
{
    return Omega<T>(a() - o.a(), b() - o.b(), c() - o.c(), d() - o.d());
}

template <typename T>
Omega<T> Omega<T>::operator*(const Omega &o) const
{
    return Omega<T>(
        a() * o.d() + b() * o.c() + c() * o.b() + d() * o.a(),
        b() * o.d() + c() * o.c() + d() * o.b() - a() * o.a(),
        c() * o.d() + d() * o.c() - a() * o.b() - b() * o.a(),
        d() * o.d() - a() * o.c() - b() * o.b() - c() * o.a());
}

template <typename T>
Omega<T> Omega<T>::operator-() const
{

    return Omega<T>(-a(), -b(), -c(), -d());
}

template <typename T>
Omega<T> Omega<T>::abs() const
{
    return copy();
}

template <typename T>
int Omega<T>::signum() const
{
    return 1;
}

template <typename T>
Omega<T> Omega<T>::adj() const
{
    return Omega(-ring::adj<T>(c()), -ring::adj<T>(b()), -ring::adj<T>(a()), ring::adj<T>(d()));
}

template <typename T>
Omega<T> Omega<T>::adj2() const
{
    return Omega(-ring::adj2<T>(a()), ring::adj2<T>(b()), -ring::adj2<T>(c()), ring::adj2<T>(d()));
}

template <typename T>
Omega<T> Omega<T>::recip() const
{
    static_assert(std::is_same<T, Real>::value || std::is_same<T, Rational>::value,
                  "recip can only be called with T = Real or T = Rational.");
    assert((*this) != 0);
    Omega<T> x1 = Omega<T>(-c(), -b(), -a(), d());
    Omega<T> x2 = Omega<T>(-a(), b(), -c(), d());
    Omega<T> x3 = Omega<T>(c(), -b(), a(), d());
    T sumSquares = a() * a() + b() * b() + c() * c() + d() * d();
    T sumProds = a() * b() + b() * c() + c() * d() - d() * a();
    T denom = sumSquares * sumSquares - 2 * sumProds * sumProds;
    return x1 * x2 * x3 * Omega(0, 0, 0, ring::recip(denom));
}

template <typename T>
Integer Omega<T>::norm() const
{
    Integer nA = ring::norm<T>(a());
    Integer nB = ring::norm<T>(b());
    Integer nC = ring::norm<T>(c());
    Integer nD = ring::norm<T>(d());
    Integer sum1 = nA * nA + nB * nB + nC * nC + nD * nD;
    Integer sum2 = nA * nB + nB * nC + nC * nD - nD * nA;
    return sum1 * sum1 - 2 * sum2 * sum2;
}

template <typename T>
QOmega Omega<T>::toQOmega() const
{
    QOmega o = ring::omega<QOmega>();
    QOmega o2 = o * o;
    QOmega o3 = o2 * o;
    return o3 * ring::toQOmega(a()) + o2 * ring::toQOmega(b()) + o * ring::toQOmega(c()) + ring::toQOmega(d());
}

template <typename T>
std::string Omega<T>::toString() const
{
    return "Omega(" + ring::toString(a()) + ", " + ring::toString(b()) + ", " + ring::toString(c()) + ", " + ring::toString(d()) + ")";
}

template <typename T>
void Omega<T>::print(std::string prefix) const
{
    std::cout << prefix << ": " << toString() << std::endl;
}

template <typename T>
Omega<T> Omega<T>::half()
{
    return Omega<T>(0, 0, 0, ring::half<T>());
}

template <typename T>
Omega<T> Omega<T>::rootTwo()
{
    return Omega<T>(-1, 0, 1, 0);
}

template <typename T>
Omega<T> Omega<T>::rootHalf()
{
    return Omega<T>(-ring::half<T>(), 0, ring::half<T>(), 0);
}

template <typename T>
Omega<T> Omega<T>::i()
{
    return Omega<T>(0, 1, 0, 0);
}

template <typename T>
Omega<T> Omega<T>::omega()
{
    return Omega<T>(0, 0, 1, 0);
}

template <typename T>
Omega<T> Omega<T>::fromInteger(int n)
{
    return Omega<T>(0, 0, 0, n);
}

template <typename T>
Omega<T> Omega<T>::fromInteger(Integer n)
{
    return Omega<T>(0, 0, 0, ring::fromInteger<T>(n));
}

template <typename T>
Omega<T> Omega<T>::fromRational(Rational r)
{
    Integer numerator = r.get_num();
    Integer denominator = r.get_den();
    return ring::fromInteger<Omega<T>>(numerator) * ring::fromInteger<Omega<T>>(denominator).recip();
}

template <>
QOmega::Omega(Rational a, Rational b, Rational c, Rational d)
{
    // Make sure Rationals are in canonical form.
    a.canonicalize();
    b.canonicalize();
    c.canonicalize();
    d.canonicalize();
    a_ = a;
    b_ = b;
    c_ = c;
    d_ = d;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Omega<T> &o)
{
    os << o.toString();
    return os;
}

std::ostream& operator<<(std::ostream& os, const Real &r)
{
    os << ring::toString(r);
    return os;
}