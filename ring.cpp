#include "ring.h"
#include <cmath>
#include <iostream>
#include <gmpxx.h>

namespace ring
{
    int mpzToInt(Integer z)
    {
        assert(z.fits_sint_p());
        return z.get_si();
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
    int sign(double a)
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

    mp_bitcnt_t lobit(Integer n)
    {
        if (n == 0)
        {
            return -1;
        }
        return mpz_scan1(n.get_mpz_t(), 0);
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
        assert(n >= 0);
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
    Rational fromInteger(int arg) { return Rational(arg); }

    template <typename T>
    T fromInteger(Integer arg)
    {
        return fromInteger<T>(mpzToInt(arg));
    }

    template <typename T>
    T recip(T arg);

    template <>
    double recip(double arg)
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
    double fromRational(Rational r)
    {
        return mpq_get_d(r.get_mpq_t());
    }

    template <typename T>
    T powNonNeg(T base, int exp)
    {
        assert(exp >= 0);
        if (exp == 0)
        {
            return 1;
        }
        T result = base;
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
    T half()
    {
        return T::half();
    }

    template <>
    double half() { return 0.5; }

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
    T fromDyadic(Dyadic<Integer> d)
    {
        std::tuple<Integer, Integer> decomposed = d.decomposeDyadic();
        Integer a = std::get<0>(decomposed);
        Integer n = std::get<1>(decomposed);
        if (n >= 0)
        {
            return fromInteger<T>(a) * powNonNeg<T>(half<T>(), mpzToInt(n));
        }
        else
        {
            return fromInteger<T>(a) * fromInteger<T>(exp2<Integer>(mpzToInt(-n)));
        }
    }

    template <typename T>
    T rootTwo()
    {
        return T::rootTwo();
    }

    template <>
    double rootTwo() { return sqrt(2); }

    template <typename T>
    T fromZRootTwo(ZRootTwo arg)
    {
        return fromInteger<T>(arg.a) + rootTwo<T>() * fromInteger<T>(arg.b);
    }

    template <typename T>
    T rootHalf()
    {
        return T::rootHalf();
    }

    template <>
    double rootHalf()
    {
        return sqrt(0.5);
    }

    template <typename T>
    T fromDRootTwo(DRootTwo arg)
    {
        return fromDyadic<T>(arg.a) + rootTwo<T>() * fromDyadic<T>(arg.b);
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
    double adj(double arg) { return arg; }

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

    Integer floor_of(Rational arg)
    {
        Integer result;
        mpz_fdiv_q(result.get_mpz_t(), arg.get_num_mpz_t(), arg.get_den_mpz_t());
        return result;
    }

    Integer floor_of(Integer arg)
    {
        return arg;
    }

    Integer floor_of(QRootTwo arg)
    {
        Integer a = floor_of(arg.a);
        Integer b = intsqrt(floor_of(2 * arg.b * arg.b));
        Integer rInt = (arg.b >= 0) ? (Integer(a + b)) : (Integer(a - b));
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
        std::optional<U> aDyadic = ring::maybeDyadic<T, U>(arg.a);
        std::optional<U> bDyadic = ring::maybeDyadic<T, U>(arg.b);
        if (aDyadic.has_value() && bDyadic.has_value())
        {
            return RootTwo<U>(aDyadic.value(), bDyadic.value());
        }
        return std::nullopt;
    }

    template <typename T, typename U>
    std::optional<Complex<U>> maybeDyadic(Complex<T> arg)
    {
        std::optional<U> aDyadic = ring::maybeDyadic<T, U>(arg.a);
        std::optional<U> bDyadic = ring::maybeDyadic<T, U>(arg.b);
        if (aDyadic.has_value() && bDyadic.has_value())
        {
            return Complex<U>(aDyadic.value(), bDyadic.value());
        }
        return std::nullopt;
    }

    template <typename T, typename U>
    std::optional<Omega<U>> maybeDyadic(Omega<T> arg)
    {
        std::optional<U> aDyadic = ring::maybeDyadic<T, U>(arg.a);
        std::optional<U> bDyadic = ring::maybeDyadic<T, U>(arg.b);
        std::optional<U> cDyadic = ring::maybeDyadic<T, U>(arg.c);
        std::optional<U> dDyadic = ring::maybeDyadic<T, U>(arg.d);
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
        return arg.a;
    }

    template <typename T>
    RootTwo<T> real(Omega<T> arg)
    {
        return RootTwo<T>(arg.d, half<T>() * (arg.c - arg.a));
    }

    template <typename T>
    QOmega toQOmega(T arg)
    {
        return arg.toQOmega();
    }

    template <>
    QOmega toQOmega<int>(int arg) { return fromInteger<Omega<Rational>>(arg); }

    template <>
    QOmega toQOmega<Integer>(Integer arg) { return fromInteger<Omega<Rational>>(arg); }

    template <>
    QOmega toQOmega<Rational>(Rational arg) { return fromRational<Omega<Rational>>(arg); }

    template <typename T>
    T fromQRootTwo(QRootTwo q)
    {
        return fromRational<T>(q.a) + rootTwo<T>() * fromRational<T>(q.b);
    }

    template <typename T>
    T fromZComplex(ZComplex z)
    {
        return fromInteger<T>(z.a) + i<T>() * fromInteger<T>(z.b);
    }

    template <typename T>
    T fromDComplex(DComplex d)
    {
        return fromDyadic<T>(d.a) + i<T>() * fromDyadic<T>(d.b);
    }

    template <typename T>
    T fromQComplex(QComplex q)
    {
        return fromRational<T>(q.a) + i<T>() * fromRational<T>(q.b);
    }

    template <typename T>
    T fromDRComplex(DRComplex d)
    {
        return fromDRootTwo<T>(d.a) + i<T>() * fromDRootTwo<T>(d.b);
    }

    template <typename T>
    T fromQRComplex(QRComplex q)
    {
        return fromQRootTwo<T>(q.a) + i<T>() * fromQRootTwo<T>(q.b);
    }

    template <typename T>
    T fromZOmega(ZOmega z)
    {
        T o = omega<T>();
        T o2 = o * o;
        T o3 = o * o * o;
        return fromInteger<T>(z.a) * o3 + fromInteger<T>(z.b) * o2 + fromInteger<T>(z.c) * o + fromInteger<T>(z.d);
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

    template <>
    std::string toString(const Rational &arg) { return arg.get_str(); }
}

template <typename T>
Dyadic<T>::Dyadic()
{
    this->a = 0;
    this->n = 0;
}

template <typename T>
Dyadic<T>::Dyadic(int arg)
{
    this->a = arg;
    this->n = 0;
}

template <typename T>
Dyadic<T>::Dyadic(T a, T n)
{
    this->a = a;
    this->n = n;
}

template <typename T>
Dyadic<T> Dyadic<T>::copy() const
{
    return Dyadic<T>(a, n);
}

template <typename T>
bool Dyadic<T>::operator==(const Dyadic &d) const
{
    T b = d.a;
    T m = d.n;
    T k = (m > n) ? m : n;
    return (a * ring::exp2<T>(k - n)) == (b * ring::exp2<T>(k - m));
}

template <typename T>
bool Dyadic<T>::operator!=(const Dyadic &d) const
{
    return !(*this == d);
}

template <typename T>
Ordering Dyadic<T>::compare(const Dyadic &d) const
{
    T b = d.a;
    T m = d.n;
    T k = (n > m) ? n : m;
    T size1 = a * ring::exp2<T>(k - n);
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
    return this->compare(d) == LT;
}

template <typename T>
bool Dyadic<T>::operator>(const Dyadic &d) const
{
    return this->compare(d) == GT;
}

template <typename T>
bool Dyadic<T>::operator<=(const Dyadic &d) const
{
    Ordering c = this->compare(d);
    return (c == LT) || (c == EQ);
}

template <typename T>
bool Dyadic<T>::operator>=(const Dyadic &d) const
{
    Ordering c = this->compare(d);
    return (c == GT) || (c == EQ);
}

template <typename T>
Dyadic<T> Dyadic<T>::operator+(const Dyadic &d) const
{
    T b = d.a;
    T m = d.n;
    if (n < m)
    {
        T c = ring::shiftL(a, m - n) + b;
        return Dyadic(c, m);
    }
    else
    {
        T d = a + ring::shiftL(b, n - m);
        return Dyadic(d, n);
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
    T b = d.a;
    T m = d.n;
    return Dyadic(a * b, m + n);
}

template <typename T>
Dyadic<T> Dyadic<T>::operator-() const
{
    return Dyadic(-a, n);
}

template <typename T>
Dyadic<T> Dyadic<T>::abs() const
{
    if (this->compare(ring::fromInteger<Dyadic<T>>(0)) != LT)
    {
        return this->copy();
    }
    return -this->copy();
}

template <typename T>
int Dyadic<T>::signum() const
{
    Ordering comp = this->compare(ring::fromInteger<Dyadic<T>>(0));
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
    return this->copy();
}

template <typename T>
Dyadic<T> Dyadic<T>::adj2() const
{
    return this->copy();
}

template <typename T>
std::tuple<T, T> Dyadic<T>::decomposeDyadic() const
{
    if (a == 0)
    {
        return std::make_tuple(0, 0);
    }
    else
    {
        int k = ring::lobit(a);
        if (n >= k)
        {
            return std::make_tuple(ring::shiftR(a, k), n - k);
        }
        else
        {
            return std::make_tuple(ring::shiftR(a, n), 0);
        }
    }
}

template <typename T>
T Dyadic<T>::integerOfDyadic(T k) const
{
    return ring::shift(a, (k - n));
}

template <typename T>
QOmega Dyadic<T>::toQOmega() const
{
    if (n >= 0)
    {
        return ring::toQOmega<T>(a) * ring::powNonNeg<QOmega>(ring::half<QOmega>(), ring::mpzToInt(n));
    }
    return ring::toQOmega<T>(a) * ring::powNonNeg<QOmega>(QOmega(2), ring::mpzToInt(-n));
}

template <typename T>
std::string Dyadic<T>::toString() const
{
    return "Dyadic(" + ring::toString(a) + ", " + ring::toString(n) + ")";
}

template <typename T>
void Dyadic<T>::print(std::string prefix) const
{
    std::cout << prefix << ": " << this->toString() << std::endl;
}

template <typename T>
Dyadic<T> Dyadic<T>::fromInteger(int n)
{
    return Dyadic(n, 0);
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
Dyadic<Integer> ring::half()
{
    return Dyadic<Integer>(1, 1);
}

template <typename T>
RootTwo<T>::RootTwo()
{
    this->a = 0;
    this->b = 0;
}

template <typename T>
RootTwo<T>::RootTwo(int arg)
{
    this->a = arg;
    this->b = 0;
}

template <typename T>
RootTwo<T>::RootTwo(T a, T b)
{
    this->a = a;
    this->b = b;
}

template <typename T>
RootTwo<T> RootTwo<T>::copy() const
{
    return RootTwo(a, b);
}

template <typename T>
bool RootTwo<T>::operator==(const RootTwo &r) const
{
    return (this->a == r.a) && (this->b == r.b);
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
    return RootTwo(this->a + r.a, this->b + r.b);
}

template <typename T>
RootTwo<T> RootTwo<T>::operator-(const RootTwo &r) const
{
    return RootTwo(this->a - r.a, this->b - r.b);
}

template <typename T>
RootTwo<T> RootTwo<T>::operator*(const RootTwo &r) const
{
    T newA = this->a * r.a + (this->b * r.b) + (this->b * r.b);
    T newB = this->a * r.b + r.a * this->b;
    return RootTwo(newA, newB);
}

template <typename T>
RootTwo<T> RootTwo<T>::operator-() const
{
    return RootTwo(-this->a, -this->b);
}

template <typename T>
RootTwo<T> RootTwo<T>::abs() const
{
    int sign = this->signum();
    if (sign != -1)
    {
        return this->copy();
    }
    return -(*this);
}

template <typename T>
int RootTwo<T>::signum() const
{
    int sa = ring::sign(a);
    int sb = ring::sign(b);
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
    else if (sa != -1 && sb != 1 && ring::sign<T>(a * a - b * b - b * b) != -1)
    {
        return 1;
    }
    else if (sa != 1 && sb != -1 && ring::sign<T>(a * a - b * b - b * b) != 1)
    {
        return 1;
    }
    return -1;
}

template <typename T>
RootTwo<T> RootTwo<T>::adj() const
{
    return RootTwo(ring::adj(a), ring::adj(b));
}

template <typename T>
RootTwo<T> RootTwo<T>::adj2() const
{
    return RootTwo(ring::adj2(a), -ring::adj2(b));
}

template <typename T>
RootTwo<T> RootTwo<T>::recip() const
{
    T k = a * a - 2 * b * b;
    return RootTwo(a * ring::recip(k), -b * ring::recip(k));
}

template <typename T>
Integer RootTwo<T>::norm() const
{
    Integer normA = ring::norm<T>(a);
    Integer normB = ring::norm<T>(b);
    return normA * normA - 2 * normB * normB;
}

template <typename T>
QOmega RootTwo<T>::toQOmega() const
{
    return ring::toQOmega<T>(a) + ring::rootTwo<QOmega>() * ring::toQOmega<T>(b);
}

template <typename T>
std::string RootTwo<T>::toString() const
{
    return "RootTwo(" + ring::toString(a) + ", " + ring::toString(b) + ")";
}

template <typename T>
void RootTwo<T>::print(std::string prefix) const
{
    std::cout << prefix << ": " << this->toString() << std::endl;
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
RootTwo<T> RootTwo<T>::fromRational(Rational r)
{
    return RootTwo<T>(ring::fromRational<T>(r), 0);
}

template <>
RootTwo<Rational>::RootTwo(Rational a, Rational b)
{
    // Make sure numerator and denominator are in canonical form.
    a.canonicalize();
    b.canonicalize();
    this->a = a;
    this->b = b;
}

template <typename T>
Complex<T>::Complex()
{
    this->a = 0;
    this->b = 0;
}

template <typename T>
Complex<T>::Complex(int arg)
{
    this->a = arg;
    this->b = 0;
}

template <typename T>
Complex<T>::Complex(T a, T b)
{
    this->a = a;
    this->b = b;
}

template <typename T>
Complex<T> Complex<T>::copy() const
{
    return Complex(a, b);
}

template <typename T>
bool Complex<T>::operator==(const Complex &c) const
{
    return (this->a == c.a) && (this->b == c.b);
}

template <typename T>
bool Complex<T>::operator!=(const Complex &c) const
{
    return !(*this == c);
}

template <typename T>
Complex<T> Complex<T>::operator+(const Complex &c) const
{
    return Complex(this->a + c.a, this->b + c.b);
}

template <typename T>
Complex<T> Complex<T>::operator-(const Complex &c) const
{
    return Complex(this->a - c.a, this->b - c.b);
}

template <typename T>
Complex<T> Complex<T>::operator*(const Complex &c) const
{
    T newA = this->a * c.a - this->b * c.b;
    T newB = this->a * c.b + c.a * this->b;
    return Complex(newA, newB);
}

template <typename T>
Complex<T> Complex<T>::operator-() const
{
    return Complex(-this->a, -this->b);
}

template <typename T>
Complex<T> Complex<T>::abs() const
{
    return this->copy();
}

template <typename T>
int Complex<T>::signum() const
{
    return 1;
}

template <typename T>
Complex<T> Complex<T>::adj() const
{
    return Complex(ring::adj(a), -ring::adj(b));
}

template <typename T>
Complex<T> Complex<T>::adj2() const
{
    return Complex(ring::adj2(a), ring::adj2(b));
}

template <typename T>
Complex<T> Complex<T>::recip() const
{
    T d = a * a + b * b;
    return Complex(a * ring::recip(d), -b * ring::recip(d));
}

template <typename T>
Integer Complex<T>::norm() const
{
    Integer normA = ring::norm<T>(a);
    Integer normB = ring::norm<T>(b);
    return normA * normA + normB * normB;
}

template <typename T>
QOmega Complex<T>::toQOmega() const
{
    return ring::toQOmega(a) + ring::i<QOmega>() * ring::toQOmega<T>(b);
}

template <typename T>
std::string Complex<T>::toString() const
{
    return "Complex(" + ring::toString(a) + ", " + ring::toString(b) + ")";
}

template <typename T>
void Complex<T>::print(std::string prefix) const
{
    std::cout << prefix << ": " << this->toString() << std::endl;
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
Complex<T> Complex<T>::fromRational(Rational r)
{
    return Complex<T>(ring::fromRational<T>(r), 0);
}

Z2::Z2()
{
    this->mod2 = 0;
}

Z2::Z2(int arg)
{
    this->mod2 = (arg % 2) != 0;
}

Z2::Z2(bool mod2)
{
    this->mod2 = mod2;
}

Z2 Z2::copy() const
{
    return Z2(this->mod2);
}

bool Z2::operator==(const Z2 &z) const
{
    return this->mod2 == z.mod2;
}

bool Z2::operator!=(const Z2 &z) const
{
    return !(*this == z);
}

Z2 Z2::operator+(const Z2 &z) const
{
    return Z2(this->mod2 != z.mod2);
}

Z2 Z2::operator-(const Z2 &z) const
{
    return (*this) + (-z);
}

Z2 Z2::operator*(const Z2 &z) const
{
    if (this->mod2 == 0)
    {
        return Z2(0);
    }
    return z.mod2;
}

Z2 Z2::operator-() const
{
    return this->copy();
}
Z2 Z2::abs() const
{
    return this->copy();
}

int Z2::signum() const
{
    return 1;
}

Z2 Z2::adj() const
{
    return this->copy();
}

Z2 Z2::adj2() const
{
    return this->copy();
}

std::string Z2::toString() const
{
    return "Z2(" + std::to_string(this->mod2) + ")";
}

void Z2::print(std::string prefix) const
{
    std::cout << prefix << ": " << this->toString() << std::endl;
}

Z2 Z2::fromInteger(int n)
{
    return Z2((bool)(n % 2));
}

template <typename T>
Omega<T>::Omega()
{
    this->a = 0;
    this->b = 0;
    this->c = 0;
    this->d = 0;
}

template <typename T>
Omega<T>::Omega(int arg)
{
    this->a = 0;
    this->b = 0;
    this->c = 0;
    this->d = arg;
}

template <typename T>
Omega<T>::Omega(T a, T b, T c, T d)
{
    this->a = a;
    this->b = b;
    this->c = c;
    this->d = d;
}

template <typename T>
Omega<T> Omega<T>::copy() const
{
    return Omega<T>(a, b, c, d);
}

template <typename T>
bool Omega<T>::operator==(const Omega &o) const
{
    return (a == o.a) && (b == o.b) && (c == o.c) && (d == o.d);
}

template <typename T>
bool Omega<T>::operator!=(const Omega &o) const
{
    return !(*this == o);
}

template <typename T>
Omega<T> Omega<T>::operator+(const Omega &o) const
{
    return Omega<T>(a + o.a, b + o.b, c + o.c, d + o.d);
}

template <typename T>
Omega<T> Omega<T>::operator-(const Omega &o) const
{
    return Omega<T>(a - o.a, b - o.b, c - o.c, d - o.d);
}

template <typename T>
Omega<T> Omega<T>::operator*(const Omega &o) const
{
    return Omega<T>(
        a * o.d + b * o.c + c * o.b + d * o.a,
        b * o.d + c * o.c + d * o.b - a * o.a,
        c * o.d + d * o.c - a * o.b - b * o.a,
        d * o.d - a * o.c - b * o.b - c * o.a);
}

template <typename T>
Omega<T> Omega<T>::operator-() const
{
    return Omega<T>(-a, -b, -c, -d);
}

template <typename T>
Omega<T> Omega<T>::abs() const
{
    return this->copy();
}

template <typename T>
int Omega<T>::signum() const
{
    return 1;
}

template <typename T>
Omega<T> Omega<T>::adj() const
{
    return Omega(-ring::adj<T>(c), -ring::adj<T>(b), -ring::adj<T>(a), ring::adj<T>(d));
}

template <typename T>
Omega<T> Omega<T>::adj2() const
{
    return Omega(-ring::adj2<T>(a), ring::adj2<T>(b), -ring::adj2<T>(c), ring::adj2<T>(d));
}

template <typename T>
Omega<T> Omega<T>::recip() const
{
    static_assert(std::is_same<T, double>::value || std::is_same<T, Rational>::value,
                  "recip can only be called with T = double or T = Rational.");
    Omega<T> x1 = Omega<T>(-c, -b, -a, d);
    Omega<T> x2 = Omega<T>(-a, b, -c, d);
    Omega<T> x3 = Omega<T>(c, -b, a, d);
    T sumSquares = a * a + b * b + c * c + d * d;
    T sumProds = a * b + b * c + c * d - d * a;
    T denom = sumSquares * sumSquares - 2 * sumProds * sumProds;
    return x1 * x2 * x3 * Omega(0, 0, 0, ring::recip(denom));
}

template <typename T>
Integer Omega<T>::norm() const
{
    Integer nA = ring::norm<T>(a);
    Integer nB = ring::norm<T>(b);
    Integer nC = ring::norm<T>(c);
    Integer nD = ring::norm<T>(d);
    Integer sum1 = nA * nA + nB * nB + nC * nC + nD * nD;
    Integer sum2 = nA * nB + nB * nC + nC * nD - nD * nA;
    return sum1 * sum1 - 2 * sum2 * sum2;
}

template <typename T>
QOmega Omega<T>::toQOmega() const
{
    QOmega o = ring::omega<QOmega>();
    QOmega o2 = o * o;
    QOmega o3 = o * o * o;
    return o3 * ring::toQOmega(a) + o2 * ring::toQOmega(b) + o * ring::toQOmega(c) + ring::toQOmega(d);
}

template <typename T>
std::string Omega<T>::toString() const
{
    return "Omega(" + ring::toString(a) + ", " + ring::toString(b) + ", " + ring::toString(c) + ", " + ring::toString(d) + ")";
}

template <typename T>
void Omega<T>::print(std::string prefix) const
{
    std::cout << prefix << ": " << this->toString() << std::endl;
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
Omega<T> Omega<T>::fromRational(Rational r)
{
    Integer numerator = r.get_num();
    Integer denominator = r.get_den();
    return ring::fromInteger<Omega<T>>(numerator) * ring::fromInteger<Omega<T>>(denominator).recip();
}
