#pragma once
#include "types.h"
#include <iostream>
#include <gmpxx.h>
#include <optional>

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
using CReal = Complex<Real>;

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
    Dyadic &operator+=(const Dyadic &d);
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
    std::tuple<T, T> decompose_dyadic() const;
    T integer_of_dyadic(T k) const;
    QOmega toQOmega() const;
    std::string to_string() const;
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
    RootTwo &operator+=(const RootTwo &r);
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
    std::string to_string() const;
    void print(std::string prefix) const;
    static RootTwo half();
    static RootTwo roottwo();
    static RootTwo roothalf();
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
    Complex &operator+=(const Complex &c);
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
    Integer denomexp() const;
    Complex denomexp_factor(Integer k) const;
    QOmega toQOmega() const;
    std::string to_string() const;
    void print(std::string prefix) const;
    static Complex half();
    static Complex roottwo();
    static Complex roothalf();
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
    std::string to_string() const;
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
    Omega &operator+=(const Omega &o);
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
    std::string to_string() const;
    void print(std::string prefix) const;
    static Omega half();
    static Omega roottwo();
    static Omega roothalf();
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

    template <typename T>
    T abs(T arg);

    template <>
    Integer abs(Integer arg);

    bool get_bit(int n, int bit);

    int lobit(int n);

    int lobit(Integer n);

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
    Real fromInteger(int arg);

    template <>
    Rational fromInteger(int arg);

    template <typename T>
    T fromInteger(Integer arg);

    template <>
    Integer fromInteger(Integer arg);

    template <>
    double fromInteger(Integer arg);

    template <>
    Real fromInteger(Integer arg);

    template <>
    Rational fromInteger(Integer arg);

    template <typename T>
    T recip(T arg);

    template <>
    Real recip(Real arg);

    template <>
    Rational recip(Rational arg);

    template <typename T>
    T fromRational(Rational r);

    template <>
    Rational fromRational(Rational r);

    template <typename T>
    T pow_non_neg(T base, int exp);

    template <typename T>
    T pow_non_neg(T base, Integer exp);

    template <typename T>
    T pow_int(T base, int exp);

    template <typename T>
    T pow_int(T base, Integer exp);

    template <typename T>
    T half();

    template <>
    Real half();

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
    T roottwo();

    template <>
    Real roottwo();

    template <typename T>
    T fromZRootTwo(ZRootTwo arg);

    template <typename T>
    T roothalf();

    template <>
    Real roothalf();

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
    Real adj(Real arg);

    template <>
    Rational adj(Rational arg);

    template <>
    Real adj(Real arg);

    template <typename T>
    T adj2(T arg);

    template <>
    int adj2(int arg);

    template <>
    Integer adj2(Integer arg);

    template <>
    Rational adj2(Rational arg);

    template <>
    Real adj2(Real arg);

    template <typename T>
    Integer norm(T arg);

    template <>
    Integer norm(int arg);

    template <>
    Integer norm(Integer arg);

    Integer floor_of(double arg);

    Integer floor_of(Real arg);

    Integer floor_of(Rational arg);

    Integer floor_of(Integer arg);

    Integer floor_of(QRootTwo arg);

    Integer ceiling_of(double arg);

    Integer ceiling_of(Real arg);

    Integer ceiling_of(Rational arg);

    Integer ceiling_of(Integer arg);

    Integer ceiling_of(QRootTwo arg);

    template <typename T>
    T omega();

    template <typename T, typename U>
    std::optional<U> maybe_dyadic(T arg);

    template <typename T, typename U>
    std::optional<RootTwo<U>> maybe_dyadic(RootTwo<T> arg);

    template <typename T, typename U>
    std::optional<Complex<U>> maybe_dyadic(Complex<T> arg);

    template <typename T, typename U>
    std::optional<Omega<U>> maybe_dyadic(Omega<T> arg);

    template <>
    std::optional<ZDyadic> maybe_dyadic(ZDyadic arg);

    template <>
    std::optional<ZDyadic> maybe_dyadic(Rational r);

    template <typename T, typename U>
    U to_dyadic(T arg);

    template <typename T, typename U>
    RootTwo<U> to_dyadic(RootTwo<T> arg);

    template <typename T, typename U>
    Complex<U> to_dyadic(Complex<T> arg);

    template <typename T, typename U>
    Omega<U> to_dyadic(Omega<T> arg);

    template <typename T>
    T real(Complex<T> arg);

    template <typename T>
    RootTwo<T> real(Omega<T> arg);

    template <typename T, typename U>
    T from_whole(U arg);

    template <>
    ZDyadic from_whole(Integer arg);

    template <>
    DOmega from_whole(ZOmega arg);

    template <>
    DRootTwo from_whole(ZRootTwo arg);

    template <typename T, typename U>
    U to_whole(T arg);

    template <>
    Integer to_whole(ZDyadic arg);

    template <>
    ZOmega to_whole(DOmega arg);

    template <>
    ZRootTwo to_whole(DRootTwo arg);

    template <typename T>
    Integer denomexp(T arg);

    template <>
    Integer denomexp(DOmega arg);

    template <>
    Integer denomexp(DRootTwo arg);

    template <typename T>
    T denomexp_factor(T arg, Integer k);

    template <>
    DOmega denomexp_factor(DOmega arg, Integer k);

    template <>
    DRootTwo denomexp_factor(DRootTwo arg, Integer k);

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

    std::optional<ZRootTwo> zroottwo_root(ZRootTwo arg);

    ZRootTwo zroottwo_of_zomega(ZOmega arg);

    template <typename T>
    std::string to_string(const T &arg);

    template <>
    std::string to_string(const int &arg);

    template <>
    std::string to_string(const Integer &arg);

    template <>
    std::string to_string(const double &arg);

    std::string to_string(const Real &arg);

    template <>
    std::string to_string(const Rational &arg);
}

#include "ring.cpp"