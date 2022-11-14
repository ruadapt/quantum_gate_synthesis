#pragma once
#include <iostream>
#include <gmpxx.h>
#include <optional>

typedef mpz_class Integer;
typedef mpq_class Rational;

template <typename T>
class Dyadic;

template <typename T>
class RootTwo;

template <typename T>
class Complex;

class Z2;

template <typename T>
class Omega;

using ZDyadic = Dyadic<Integer>;

using ZRootTwo = RootTwo<Integer>;
using DRootTwo = RootTwo<ZDyadic>;
using QRootTwo = RootTwo<Rational>;

using ZComplex = Complex<Integer>;
using DComplex = Complex<ZDyadic>;
using QComplex = Complex<Rational>;
using DRComplex = Complex<DRootTwo>;
using QRComplex = Complex<QRootTwo>;
using CDouble = Complex<double>;

using ZOmega = Omega<Integer>;
using DOmega = Omega<ZDyadic>;
using QOmega = Omega<Rational>;

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
    Dyadic();
    Dyadic(int arg);
    Dyadic(Integer arg);
    Dyadic(T a, T n);
    T a() const;
    T n() const;
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
    QOmega toQOmega() const;
    std::string toString() const;
    void print(std::string prefix) const;
    static Dyadic fromInteger(int n);
    static Dyadic fromInteger(Integer n);
    static Dyadic fromDyadic(const Dyadic &d);
    static Dyadic half();

private:
    T a_;
    T n_;
};

template <typename T>
class RootTwo
{
public:
    RootTwo();
    RootTwo(int arg);
    RootTwo(Integer arg);
    RootTwo(T a, T b);
    T a() const;
    T b() const;
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
    QOmega toQOmega() const;
    std::string toString() const;
    void print(std::string prefix) const;
    static RootTwo half();
    static RootTwo rootTwo();
    static RootTwo rootHalf();
    static RootTwo i();
    static RootTwo omega();
    static RootTwo fromInteger(int n);
    static RootTwo fromInteger(Integer n);
    static RootTwo fromRational(Rational r);

private:
    T a_;
    T b_;
};

template <>
RootTwo<Rational>::RootTwo(Rational a, Rational b);

template <typename T>
class Complex
{
public:
    Complex();
    Complex(int arg);
    Complex(Integer arg);
    Complex(T a, T b);
    T a() const;
    T b() const;
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
    Integer denomExp() const;
    Complex denomExpFactor(Integer k) const;
    QOmega toQOmega() const;
    std::string toString() const;
    void print(std::string prefix) const;
    static Complex half();
    static Complex rootTwo();
    static Complex rootHalf();
    static Complex i();
    static Complex omega();
    static Complex fromInteger(int n);
    static Complex fromInteger(Integer n);
    static Complex fromRational(Rational r);

private:
    T a_;
    T b_;
};

template <>
QComplex::Complex(Rational a, Rational b);

class Z2
{
public:
    Z2();
    Z2(int arg);
    Z2(Integer arg);
    Z2(bool mod2);
    bool mod2() const;
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
    static Z2 fromInteger(Integer n);

private:
    bool mod2_;
};

template <typename T>
class Omega
{
public:
    Omega();
    Omega(int arg);
    Omega(Integer arg);
    Omega(T a, T b, T c, T d);
    T a() const;
    T b() const;
    T c() const;
    T d() const;
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
    QOmega toQOmega() const;
    std::string toString() const;
    void print(std::string prefix) const;
    static Omega half();
    static Omega rootTwo();
    static Omega rootHalf();
    static Omega i();
    static Omega omega();
    static Omega fromInteger(int n);
    static Omega fromInteger(Integer n);
    static Omega fromRational(Rational r);

private:
    T a_;
    T b_;
    T c_;
    T d_;
};

template <>
QOmega::Omega(Rational a, Rational b, Rational c, Rational d);

namespace ring
{
    bool even(int n);

    bool even(Integer n);

    bool odd(int n);

    bool odd(Integer n);

    int div(int a, int b);

    Integer div(Integer a, Integer b);

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

    std::optional<Integer> log2(Integer n);

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

    template <>
    Integer fromInteger(Integer arg);

    template <>
    double fromInteger(Integer arg);

    template <>
    Rational fromInteger(Integer arg);

    template <typename T>
    T recip(T arg);

    template <>
    double recip(double arg);

    template <>
    Rational recip(Rational arg);

    template <typename T>
    T fromRational(Rational r);

    template <>
    Rational fromRational(Rational r);

    template <typename T>
    T powNonNeg(T base, int exp);

    template <typename T>
    T powInt(T base, int exp);

    template <typename T>
    T half();

    template <>
    double half();

    template <>
    Rational half();

    template <>
    Dyadic<int> half();

    template <>
    ZDyadic half();

    template <typename T>
    T fromDyadic(Dyadic<int> d);

    template <typename T>
    T fromDyadic(ZDyadic d);

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
    T omega();

    template <typename T, typename U>
    std::optional<U> maybeDyadic(T arg);

    template <typename T, typename U>
    std::optional<RootTwo<U>> maybeDyadic(RootTwo<T> arg);

    template <typename T, typename U>
    std::optional<Complex<U>> maybeDyadic(Complex<T> arg);

    template <typename T, typename U>
    std::optional<Omega<U>> maybeDyadic(Omega<T> arg);

    template <>
    std::optional<ZDyadic> maybeDyadic(ZDyadic arg);

    template <>
    std::optional<ZDyadic> maybeDyadic(Rational r);

    template <typename T, typename U>
    U toDyadic(T arg);

    template <typename T, typename U>
    RootTwo<U> toDyadic(RootTwo<T> arg);

    template <typename T, typename U>
    Complex<U> toDyadic(Complex<T> arg);

    template <typename T, typename U>
    Omega<U> toDyadic(Omega<T> arg);

    template <typename T>
    T real(Complex<T> arg);

    template <typename T>
    RootTwo<T> real(Omega<T> arg);

    template <typename T, typename U>
    T fromWhole(U arg);

    template <>
    ZDyadic fromWhole(Integer arg);

    template <>
    DOmega fromWhole(ZOmega arg);

    template <>
    DRootTwo fromWhole(ZRootTwo arg);

    template <typename T, typename U>
    U toWhole(T arg);

    template <>
    Integer toWhole(ZDyadic arg);

    template <>
    ZOmega toWhole(DOmega arg);

    template <>
    ZRootTwo toWhole(DRootTwo arg);

    template <typename T>
    Integer denomExp(T arg);

    template <>
    Integer denomExp(DOmega arg);

    template <>
    Integer denomExp(DRootTwo arg);

    template <typename T>
    T denomExpFactor(T arg, Integer k);

    template <>
    DOmega denomExpFactor(DOmega arg, Integer k);

    template <>
    DRootTwo denomExpFactor(DRootTwo arg, Integer k);

    template <typename T>
    QOmega toQOmega(T arg);

    template <>
    QOmega toQOmega<int>(int arg);

    template <>
    QOmega toQOmega(Integer arg);

    template <>
    QOmega toQOmega(Rational arg);

    Z2 parity(int arg);

    Z2 parity(Integer arg);

    Z2 parity(ZRootTwo arg);

    template <typename T>
    T fromQRootTwo(QRootTwo q);

    template <typename T>
    T fromZComplex(ZComplex z);

    template <typename T>
    T fromDComplex(DComplex d);

    template <typename T>
    T fromQComplex(QComplex q);

    template <typename T>
    T fromDRComplex(DRComplex d);

    template <typename T>
    T fromQRComplex(QRComplex q);

    template <typename T>
    T fromZOmega(ZOmega z);

    template <typename T>
    T fromDOmega(DOmega z);

    template <typename T>
    T fromQOmega(QOmega z);

    std::optional<ZRootTwo> zRootTwoRoot(ZRootTwo arg);

    ZRootTwo zRootTwoOfZOmega(ZOmega arg);

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