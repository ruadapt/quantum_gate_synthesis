#include "ring.h"
#include <iostream>
#include <gmpxx.h>
namespace ring
{
    signed long int mpzToLongInt(Integer z)
    {
        assert(z.fits_slong_p());
        return z.get_si();
    }

    template <typename Integral>
    Integral shift(Integral x, long int bits)
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
        return shift(x, mpzToLongInt(bits));
    }

    template <>
    Integer shift(Integer x, long int bits)
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
        return shift(x, mpzToLongInt(bits));
    }

    template <typename Integral>
    Integral shiftL(Integral x, long int bits)
    {
        assert(bits >= 0);
        return shift(x, bits);
    }

    template <typename Integral>
    Integral shiftL(Integral x, Integer bits)
    {
        return shiftL(x, mpzToLongInt(bits));
    }

    template <typename Integral>
    Integral shiftR(Integral x, long int bits)
    {
        assert(bits >= 0);
        return shift(x, -bits);
    }

    template <typename Integral>
    Integral shiftR(Integral x, Integer bits)
    {
        return shiftR(x, mpzToLongInt(bits));
    }

    template <typename Integral>
    Integral exp2(long int pow)
    {
        return shiftL(1, pow);
    }

    template <typename Integral>
    Integral exp2(Integer pow)
    {
        return exp2<Integral>(mpzToLongInt(pow));
    }

    template <typename T>
    int sign(T a)
    {
        return a.signum();
    }

    template <>
    int sign(long int a)
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
    int sign(Rational a)
    {
        return sgn(a);
    }

    template <>
    int sign(int a)
    {
        return sign<long int>(a);
    }

    bool getBit(long int n, long int bit)
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
    T powNonNeg(T base, int exp)
    {
        assert(exp >= 0);
        if (exp == 0)
        {
            return T(1);
        }
        T result = base;
        if (exp == 1)
        {
            return base;
        }
        if (exp % 2 == 0)
        {
            return powNonNeg(base * base, exp / 2);
        }
        return base * powNonNeg(base * base, (exp - 1) / 2);
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
            return fromInteger<T>(a) * fromInteger<T>(exp2<T>(-n));
        }
    }

    // template <typename T>
    // T fromDyadic(Dyadic<Integer> d);

    // template <>
    // double fromDyadic(Dyadic<int> d);

    // template <>
    // double fromDyadic(Dyadic<Integer> d);

    // template <>
    // Rational fromDyadic(Dyadic<int> d);

    // template <>
    // Rational fromDyadic(Dyadic<Integer> d);

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
Dyadic<T>::Dyadic(T arg)
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
    if (this->compare(fromInteger(0)) != LT)
    {
        return this->copy();
    }
    return -this->copy();
}

template <typename T>
int Dyadic<T>::signum() const
{
    Ordering comp = this->compare(fromInteger(0));
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
    return Dyadic(T(n), T(0));
}

template <typename T>
Dyadic<T> Dyadic<T>::fromDyadic(const Dyadic &d)
{
    return d.copy();
}

template <typename T>
Dyadic<T> Dyadic<T>::half()
{
    return Dyadic(T(1), T(1));
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
    this->a = T(0);
    this->b = T(0);
}

template <typename T>
RootTwo<T>::RootTwo(T arg)
{
    this->a = arg;
    this->b = T(0);
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
    T k = pow(a, 2) - 2 * pow(b, 2);
    return RootTwo(a / k, -b / k);
}

template <typename T>
RootTwo<T> RootTwo<T>::norm() const
{
    return pow(abs(a), 2) - 2 * pow(abs(b), 2);
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
    return RootTwo<T>(ring::half<T>(), T(0));
}

template <typename T>
RootTwo<T> RootTwo<T>::rootTwo()
{
    return RootTwo<T>(T(0), T(1));
}

template <typename T>
RootTwo<T> RootTwo<T>::rootHalf()
{
    return RootTwo<T>(T(0), ring::half<T>());
}

template <typename T>
RootTwo<T> RootTwo<T>::fromInteger(int n)
{
    return RootTwo<T>(T(n), T(0));
}

template <typename T>
RootTwo<T> RootTwo<T>::fromRational(double x)
{
    return RootTwo<T>(T(x), T(0));
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

using ZRootTwo = RootTwo<Integer>;
using QRootTwo = RootTwo<Rational>;

template <typename T>
Complex<T>::Complex()
{
    this->a = T(0);
    this->b = T(0);
}

template <typename T>
Complex<T>::Complex(T arg)
{
    this->a = arg;
    this->b = T(0);
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
    return Complex(a / d, -b / d);
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
    return Complex<T>(ring::half<T>(), T(0));
}

template <typename T>
Complex<T> Complex<T>::i()
{
    return Complex<T>(T(0), T(1));
}

Z2::Z2()
{
    this->mod2 = 0;
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