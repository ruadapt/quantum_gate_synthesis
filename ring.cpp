#include <iostream>
#include <gmpxx.h>

namespace ring
{
    signed long int mpzToLongInt(mpz_class z)
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
    Integral shift(Integral x, mpz_class bits)
    {
        return shift(x, mpzToLongInt(bits));
    }

    template <>
    mpz_class shift(mpz_class x, long int bits)
    {
        mpz_class result;
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
    mpz_class shift(mpz_class x, mpz_class bits)
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
    Integral shiftL(Integral x, mpz_class bits)
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
    Integral shiftR(Integral x, mpz_class bits)
    {
        return shiftR(x, mpzToLongInt(bits));
    }

    template <typename Integral>
    Integral exp2(long int pow)
    {
        return shiftL(1, pow);
    }

    template <typename Integral>
    Integral exp2(mpz_class pow)
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
    int sign(mpz_class a)
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
    int sign(mpq_class a)
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

    mp_bitcnt_t lobit(mpz_class n)
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

    size_t hibit(mpz_class n)
    {
        // We need the - 1 because mpz_sizeinbase returns the position of the first 1
        // bit counting from 1 instead of 0.
        return mpz_sizeinbase(n.get_mpz_t(), 2) - 1;
    }

    template <typename T>
    T half();

    template <>
    double half()
    {
        return 0.5;
    }

    template <>
    mpq_class half()
    {
        return mpq_class(1, 2);
    }

    template <typename T>
    std::string toString(const T &arg)
    {
        return std::to_string(arg);
    }
}

enum Ordering
{
    LT,
    EQ,
    GT
};

template <typename T>
T subtract(T arg1, T arg2)
{
    return arg1 + (-arg2);
}

template <typename T>
class Dyadic
{
public:
    T a;
    T n;
    Dyadic(T arg)
    {
        this->a = arg;
        this->n = 0;
    }
    Dyadic()
    {
        this->a = 0;
        this->n = 0;
    }
    Dyadic(T a, T n)
    {
        this->a = a;
        this->n = n;
    }
    Dyadic copy() const
    {
        return Dyadic(a, n);
    }
    bool operator==(const Dyadic &d) const
    {
        T b = d.a;
        T m = d.n;
        T k = (m > n) ? m : n;
        return (a * ring::exp2<T>(k - n)) == (b * ring::exp2<T>(k - m));
    }
    bool operator!=(const Dyadic &d) const
    {
        return !(*this == d);
    }
    Ordering compare(const Dyadic &d) const
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
    bool operator<(const Dyadic &d) const
    {
        return this->compare(d) == LT;
    }
    bool operator>(const Dyadic &d) const
    {
        return this->compare(d) == GT;
    }
    bool operator<=(const Dyadic &d) const
    {
        Ordering c = this->compare(d);
        return (c == LT) || (c == EQ);
    }
    bool operator>=(const Dyadic &d) const
    {
        Ordering c = this->compare(d);
        return (c == GT) || (c == EQ);
    }
    Dyadic operator+(const Dyadic &d) const
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
    Dyadic operator-(const Dyadic &d) const
    {
        return (*this) + (-d);
    }
    Dyadic operator*(const Dyadic &d) const
    {
        T b = d.a;
        T m = d.n;
        return Dyadic(a * b, m + n);
    }
    Dyadic operator-() const
    {
        return Dyadic(-a, n);
    }
    Dyadic abs() const
    {
        if (this->compare(fromInteger(0)) != LT)
        {
            return this->copy();
        }
        return -this->copy();
    }
    int signum() const
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
    Dyadic adj() const
    {
        return this->copy();
    }
    Dyadic adj2() const
    {
        return this->copy();
    }
    std::tuple<T, T> decomposeDyadic() const
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
    T integerOfDyadic(T k) const
    {
        return ring::shift(a, (k - n));
    };
    std::string toString() const
    {
        return "Dyadic(" + ring::toString(a) + ", " + ring::toString(n) + ")";
    }
    void print(std::string prefix) const
    {
        std::cout << prefix << ": " << this->toString() << std::endl;
    }
    static Dyadic fromInteger(int n)
    {
        return Dyadic(T(n), T(0));
    }
    static Dyadic fromDyadic(const Dyadic &d)
    {
        return d.copy();
    }
    static Dyadic half()
    {
        return Dyadic(T(1), T(1));
    }
};

template <>
Dyadic<int> ring::half()
{
    return Dyadic<int>(1, 1);
}

template <>
Dyadic<mpz_class> ring::half()
{
    return Dyadic<mpz_class>(1, 1);
}

template <>
std::string Dyadic<mpz_class>::toString() const
{
    return "Dyadic(" + a.get_str() + ", " + n.get_str() + ")";
}

template <typename T>
class RootTwo
{
public:
    T a;
    T b;
    RootTwo(T arg)
    {
        this->a = arg;
        this->b = T(0);
    }
    RootTwo(T a, T b)
    {
        this->a = a;
        this->b = b;
    }
    RootTwo copy() const
    {
        return RootTwo(a, b);
    }
    bool operator==(const RootTwo &r) const
    {
        return (this->a == r.a) && (this->b == r.b);
    }
    bool operator!=(const RootTwo &r) const
    {
        return !(*this == r);
    }
    bool operator<=(const RootTwo &r) const
    {
        return (r - *this).signum() != -1;
    }
    bool operator<(const RootTwo &r) const
    {
        return (*this <= r) && !(*this == r);
    }
    bool operator>=(const RootTwo &r) const
    {
        return !(*this < r);
    }
    bool operator>(const RootTwo &r) const
    {
        return !(*this <= r);
    }
    RootTwo operator+(const RootTwo &r) const
    {
        return RootTwo(this->a + r.a, this->b + r.b);
    }
    RootTwo operator-(const RootTwo &r) const
    {
        return RootTwo(this->a - r.a, this->b - r.b);
    }
    RootTwo operator*(const RootTwo &r) const
    {
        T newA = this->a * r.a + (this->b * r.b) + (this->b * r.b);
        T newB = this->a * r.b + r.a * this->b;
        return RootTwo(newA, newB);
    }
    RootTwo operator-() const
    {
        return RootTwo(-this->a, -this->b);
    }
    RootTwo abs() const
    {
        int sign = this->signum();
        if (sign != -1)
        {
            return this->copy();
        }
        return -(*this);
    }
    int signum() const
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
    RootTwo recip() const
    {
        T k = pow(a, 2) - 2 * pow(b, 2);
        return RootTwo(a / k, -b / k);
    }
    RootTwo norm() const
    {
        return pow(abs(a), 2) - 2 * pow(abs(b), 2);
    }
    std::string toString() const
    {
        return "RootTwo(" + ring::toString(a) + ", " + ring::toString(b) + ")";
    }
    void print(std::string prefix) const
    {
        std::cout << prefix << ": " << this->toString() << std::endl;
    }
    static RootTwo half()
    {
        return RootTwo<T>(ring::half<T>(), T(0));
    }
    static RootTwo rootTwo()
    {
        return RootTwo<T>(T(0), T(1));
    }
    static RootTwo rootHalf()
    {
        return RootTwo<T>(T(0), ring::half<T>());
    }
    static RootTwo fromInteger(int n)
    {
        return RootTwo<T>(T(n), T(0));
    }
    static RootTwo fromRational(double x)
    {
        return RootTwo<T>(T(x), T(0));
    }
};

template <>
RootTwo<mpq_class>::RootTwo(mpq_class a, mpq_class b)
{
    // Make sure numerator and denominator are in canonical form.
    a.canonicalize();
    b.canonicalize();
    this->a = a;
    this->b = b;
}

template <>
std::string RootTwo<mpz_class>::toString() const
{
    return "RootTwo(" + a.get_str() + ", " + b.get_str() + ")";
}

template <>
std::string RootTwo<mpq_class>::toString() const
{
    return "RootTwo(" + a.get_str() + ", " + b.get_str() + ")";
}

using ZRootTwo = RootTwo<mpz_class>;
using QRootTwo = RootTwo<mpq_class>;

template <typename T>
class Complex
{
public:
    T a;
    T b;
    Complex(T arg)
    {
        this->a = arg;
        this->b = T(0);
    }
    Complex(T a, T b)
    {
        this->a = a;
        this->b = b;
    }
    Complex copy() const
    {
        return Complex(a, b);
    }
    bool operator==(const Complex &c) const
    {
        return (this->a == c.a) && (this->b == c.b);
    }
    bool operator!=(const Complex &c) const
    {
        return !(*this == c);
    }
    Complex operator+(const Complex &c) const
    {
        return Complex(this->a + c.a, this->b + c.b);
    }
    Complex operator-(const Complex &c) const
    {
        return Complex(this->a - c.a, this->b - c.b);
    }
    Complex operator*(const Complex &c) const
    {
        T newA = this->a * c.a - this->b * c.b;
        T newB = this->a * c.b + c.a * this->b;
        return Complex(newA, newB);
    }
    Complex operator-() const
    {
        return Complex(-this->a, -this->b);
    }
    Complex abs() const
    {
        return this->copy();
    }
    int signum() const
    {
        return 1;
    }
    std::string toString() const
    {
        return "Complex(" + ring::toString(a) + ", " + ring::toString(b) + ")";
    }
    void print(std::string prefix) const
    {
        std::cout << prefix << ": " << this->toString() << std::endl;
    }
    static Complex half()
    {
        return Complex<T>(ring::half<T>(), T(0));
    }
};

template <>
std::string Complex<mpz_class>::toString() const
{
    return "Complex(" + a.get_str() + ", " + b.get_str() + ")";
}

class Z2
{
public:
    bool mod2;
    Z2(bool mod2)
    {
        this->mod2 = mod2;
    }
    Z2 copy() const
    {
        return Z2(this->mod2);
    }
    bool operator==(const Z2 &z) const
    {
        return this->mod2 == z.mod2;
    }
    bool operator!=(const Z2 &z) const
    {
        return !(*this == z);
    }
    Z2 operator+(const Z2 &z) const
    {
        return Z2(this->mod2 != z.mod2);
    }
    Z2 operator-(const Z2 &z) const
    {
        return (*this) + (-z);
    }
    Z2 operator*(const Z2 &z) const
    {
        if (this->mod2 == 0)
        {
            return Z2(0);
        }
        return z.mod2;
    }
    Z2 operator-() const
    {
        return this->copy();
    }
    Z2 abs() const
    {
        return this->copy();
    }
    int signum() const
    {
        return 1;
    }
    Z2 adj() const
    {
        return this->copy();
    }
    Z2 adj2() const
    {
        return this->copy();
    }
    std::string toString() const
    {
        return "Z2(" + std::to_string(this->mod2) + ")";
    }
    void print(std::string prefix) const
    {
        std::cout << prefix << ": " << this->toString() << std::endl;
    }
    static Z2 fromInteger(int n)
    {
        return Z2((bool)(n % 2));
    }
};

template <>
std::string ring::toString(const RootTwo<int> &arg)
{
    return arg.toString();
}

template <>
std::string ring::toString(const RootTwo<long int> &arg)
{
    return arg.toString();
}

template <>
std::string ring::toString(const RootTwo<mpz_class> &arg)
{
    return arg.toString();
}

template <>
std::string ring::toString(const Dyadic<int> &arg)
{
    return arg.toString();
}

template <>
std::string ring::toString(const Dyadic<long int> &arg)
{
    return arg.toString();
}

template <>
std::string ring::toString(const Dyadic<mpz_class> &arg)
{
    return arg.toString();
}

template <>
std::string ring::toString(const Z2 &arg)
{
    return arg.toString();
}