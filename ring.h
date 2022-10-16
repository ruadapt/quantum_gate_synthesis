#pragma once
#include <iostream>
#include <gmpxx.h>

typedef mpz_class Integer;
typedef mpq_class Rational;

namespace ring
{
    signed long int mpzToLongInt(Integer z);

    template <typename Integral>
    Integral shift(Integral x, long int bits);

    template <typename Integral>
    Integral shift(Integral x, Integer bits);

    template <>
    Integer shift(Integer x, long int bits);

    template <>
    Integer shift(Integer x, Integer bits);

    template <typename Integral>
    Integral shiftL(Integral x, long int bits);

    template <typename Integral>
    Integral shiftL(Integral x, Integer bits);

    template <typename Integral>
    Integral shiftR(Integral x, long int bits);

    template <typename Integral>
    Integral shiftR(Integral x, Integer bits);

    template <typename Integral>
    Integral exp2(long int pow);

    template <typename Integral>
    Integral exp2(Integer pow);

    template <typename T>
    int sign(T a);

    template <>
    int sign(long int a);

    template <>
    int sign(Integer a);

    template <>
    int sign(Rational a);

    template <>
    int sign(int a);

    bool getBit(long int n, long int bit);

    int lobit(int n);

    mp_bitcnt_t lobit(Integer n);

    int hibit(int n);

    size_t hibit(Integer n);

    template <typename T>
    T half();

    template <>
    double half();

    template <>
    Rational half();

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

template <>
Dyadic<int> ring::half();

template <>
Dyadic<Integer> ring::half();

template <typename T>
class RootTwo
{
public:
    T a;
    T b;
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
    RootTwo norm() const;
    std::string toString() const;
    void print(std::string prefix) const;
    static RootTwo half();
    static RootTwo rootTwo();
    static RootTwo rootHalf();
    static RootTwo fromInteger(int n);
    static RootTwo fromRational(double x);
};

template <>
RootTwo<Rational>::RootTwo(Rational a, Rational b);

using ZRootTwo = RootTwo<Integer>;
using QRootTwo = RootTwo<Rational>;

template <typename T>
class Complex
{
public:
    T a;
    T b;
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
    std::string toString() const;
    void print(std::string prefix) const;
    static Complex half();
};

class Z2
{
public:
    bool mod2;
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

#include "ring.cpp"