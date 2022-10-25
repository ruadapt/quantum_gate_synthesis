#pragma once
#include <iostream>
#include <gmpxx.h>

typedef mpz_class Integer;
typedef mpq_class Rational;

enum Ordering
{
    LT,
    EQ,
    GT
};

template <typename T>
class Dyadic
{
public:
    T a;
    T n;
    Dyadic();
    Dyadic(T arg);
    Dyadic(T a, T n);
    Dyadic copy() const;
    bool operator==(const Dyadic &d) const;
    bool operator!=(const Dyadic &d) const;
    Ordering compare(const Dyadic &d) const;
    bool operator<(const Dyadic &d) const;
    bool operator>(const Dyadic &d) const;
    bool operator<=(const Dyadic &d) const;
    bool operator>=(const Dyadic &d) const;
    Dyadic operator+(const Dyadic &d) const;
    Dyadic operator-(const Dyadic &d) const;
    Dyadic operator*(const Dyadic &d) const;
    Dyadic operator-() const;
    Dyadic abs() const;
    int signum() const;
    Dyadic adj() const;
    Dyadic adj2() const;
    std::tuple<T, T> decomposeDyadic() const;
    T integerOfDyadic(T k) const;
    std::string toString() const;
    void print(std::string prefix) const;
    static Dyadic fromInteger(int n);
    static Dyadic fromDyadic(const Dyadic &d);
    static Dyadic half();
};

using ZDyadic = Dyadic<Integer>;

template <typename T>
class RootTwo
{
public:
    T a;
    T b;
    RootTwo();
    RootTwo(T arg);
    RootTwo(T a, T b);
    RootTwo copy() const;
    bool operator==(const RootTwo &r) const;
    bool operator!=(const RootTwo &r) const;
    bool operator<=(const RootTwo &r) const;
    bool operator<(const RootTwo &r) const;
    bool operator>=(const RootTwo &r) const;
    bool operator>(const RootTwo &r) const;
    RootTwo operator+(const RootTwo &r) const;
    RootTwo operator-(const RootTwo &r) const;
    RootTwo operator*(const RootTwo &r) const;
    RootTwo operator-() const;
    RootTwo abs() const;
    int signum() const;
    RootTwo adj() const;
    RootTwo adj2() const;
    RootTwo recip() const;
    Integer norm() const;
    std::string toString() const;
    void print(std::string prefix) const;
    static RootTwo half();
    static RootTwo rootTwo();
    static RootTwo rootHalf();
    static RootTwo i();
    static RootTwo fromInteger(int n);
    static RootTwo fromRational(double x);
};

template <>
RootTwo<Rational>::RootTwo(Rational a, Rational b);

using ZRootTwo = RootTwo<Integer>;
using DRootTwo = RootTwo<ZDyadic>;
using QRootTwo = RootTwo<Rational>;

template <typename T>
class Complex
{
public:
    T a;
    T b;
    Complex();
    Complex(T arg);
    Complex(T a, T b);
    Complex copy() const;
    bool operator==(const Complex &c) const;
    bool operator!=(const Complex &c) const;
    Complex operator+(const Complex &c) const;
    Complex operator-(const Complex &c) const;
    Complex operator*(const Complex &c) const;
    Complex operator-() const;
    Complex abs() const;
    int signum() const;
    Complex adj() const;
    Complex adj2() const;
    Complex recip() const;
    Integer norm() const;
    std::string toString() const;
    void print(std::string prefix) const;
    static Complex half();
    static Complex rootTwo();
    static Complex rootHalf();
    static Complex i();
    static Complex fromInteger(int n);
};
class Z2
{
public:
    bool mod2;
    Z2();
    Z2(bool mod2);
    Z2 copy() const;
    bool operator==(const Z2 &z) const;
    bool operator!=(const Z2 &z) const;
    Z2 operator+(const Z2 &z) const;
    Z2 operator-(const Z2 &z) const;
    Z2 operator*(const Z2 &z) const;
    Z2 operator-() const;
    Z2 abs() const;
    int signum() const;
    Z2 adj() const;
    Z2 adj2() const;
    std::string toString() const;
    void print(std::string prefix) const;
    static Z2 fromInteger(int n);
};

using ZComplex = Complex<Integer>;
using DComplex = Complex<Dyadic<Integer>>;
using QComplex = Complex<Rational>;
using DRComplex = Complex<DRootTwo>;
using QRComplex = Complex<QRootTwo>;
using CDouble = Complex<double>;

template <typename T>
class Omega
{
public:
    T a;
    T b;
    T c;
    T d;
    Omega();
    Omega(T a, T b, T c, T d);
    Omega copy() const;
    bool operator==(const Omega &o) const;
    bool operator!=(const Omega &o) const;
    Omega operator+(const Omega &o) const;
    Omega operator-(const Omega &o) const;
    Omega operator*(const Omega &o) const;
    Omega operator-() const;
    Omega abs() const;
    int signum() const;
    Omega adj() const;
    Omega adj2() const;
    Omega recip() const;
    Integer norm() const;
    std::string toString() const;
    void print(std::string prefix) const;
};

namespace ring
{
    int mpzToInt(Integer z);

    template <typename Integral>
    Integral shift(Integral x, int bits);

    template <typename Integral>
    Integral shift(Integral x, Integer bits);

    template <>
    Integer shift(Integer x, int bits);

    template <>
    Integer shift(Integer x, Integer bits);

    template <typename Integral>
    Integral shiftL(Integral x, int bits);

    template <typename Integral>
    Integral shiftL(Integral x, Integer bits);

    template <typename Integral>
    Integral shiftR(Integral x, int bits);

    template <typename Integral>
    Integral shiftR(Integral x, Integer bits);

    template <typename Integral>
    Integral exp2(int pow);

    template <typename Integral>
    Integral exp2(Integer pow);

    template <typename T>
    int sign(T a);

    template <>
    int sign(int a);

    template <>
    int sign(Integer a);

    template <>
    int sign(Rational a);

    template <>
    int sign(int a);

    bool getBit(int n, int bit);

    int lobit(int n);

    mp_bitcnt_t lobit(Integer n);

    int hibit(int n);

    size_t hibit(Integer n);

    Integer intsqrt(Integer n);

    template <typename T>
    T fromInteger(int arg);

    template <>
    int fromInteger(int arg);

    template <>
    Integer fromInteger(int arg);

    template <>
    double fromInteger(int arg);

    template <>
    Rational fromInteger(int arg);

    template <typename T>
    T fromInteger(Integer arg);

    template <typename T>
    T powNonNeg(T base, int exp);

    template <typename T>
    T half();

    template <>
    double half();

    template <>
    Rational half();

    template <>
    Dyadic<int> half();

    template <>
    Dyadic<Integer> half();

    template <typename T>
    T fromDyadic(Dyadic<int> d);

    template <typename T>
    T fromDyadic(Dyadic<Integer> d);

    template <typename T>
    T rootTwo();

    template <>
    double rootTwo();

    template <typename T>
    T fromZRootTwo(ZRootTwo arg);

    template <typename T>
    T rootHalf();

    template <>
    double rootHalf();

    template <typename T>
    T fromDRootTwo(DRootTwo arg);

    template <typename T>
    T i();

    template <typename T>
    T adj(T arg);

    template <>
    int adj(int arg);

    template <>
    Integer adj(Integer arg);

    template <>
    double adj(double arg);

    template <>
    Rational adj(Rational arg);

    template <typename T>
    T adj2(T arg);

    template <>
    int adj2(int arg);

    template <>
    Integer adj2(Integer arg);

    template <>
    Rational adj2(Rational arg);

    template <typename T>
    Integer norm(T arg);

    template <>
    Integer norm(int arg);

    template <>
    Integer norm(Integer arg);

    Integer floor_of(double arg);

    Integer floor_of(Rational arg);

    Integer floor_of(Integer arg);

    Integer floor_of(QRootTwo arg);

    Integer ceiling_of(double arg);

    Integer ceiling_of(Rational arg);

    Integer ceiling_of(Integer arg);

    Integer ceiling_of(QRootTwo arg);

    template <typename T>
    std::string toString(const T &arg);

    template <>
    std::string toString(const int &arg);

    template <>
    std::string toString(const Integer &arg);

    template <>
    std::string toString(const double &arg);

    template <>
    std::string toString(const Rational &arg);
}

#include "ring.cpp"