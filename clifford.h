#pragma once
#include "types.h"

enum Axis
{
    Axis_I,
    Axis_H,
    Axis_SH
};

class Clifford
{
public:
    Clifford(): a_{0}, b_{0}, c_{0}, d_{0}
    {
    }
    Clifford(int a, int b, int c, int d) : a_{a}, b_{b}, c_{c}, d_{d}
    {
    }
    bool operator==(const Clifford &c) const
    {
        return (a_ == c.a_) && (b_ == c.b_) && (c_ == c.c_) && (d_ == c.d_);
    }
    int a() const { return a_; }
    int b() const { return b_; }
    int c() const { return c_; }
    int d() const { return d_; }
    std::string to_string() const
    {
        return "C" + std::to_string(a()) + std::to_string(b()) + std::to_string(c()) + std::to_string(d());
    }

private:
    int a_;
    int b_;
    int c_;
    int d_;
};

std::ostream &operator<<(std::ostream &os, const Clifford &c)
{
    os << c.to_string();
    return os;
}

namespace clifford
{
    Clifford clifford_X();

    Clifford clifford_Y();

    Clifford clifford_Z();

    Clifford clifford_H();

    Clifford clifford_S();

    Clifford clifford_SH();

    Clifford clifford_E();

    Clifford clifford_W();

    template <typename T>
    Clifford to_clifford(T arg);

    template <>
    Clifford to_clifford(Clifford c);

    template <>
    Clifford to_clifford(char c);

    template <>
    Clifford to_clifford(std::string s);

    template <>
    Clifford to_clifford(Axis arg);

    template <typename T>
    Tup4<int> clifford_decompose(T m);

    template <typename T>
    std::tuple<Axis, int, int, int> clifford_decompose_coset(T u);

    Clifford clifford_id();

    Clifford clifford_mult(Clifford u1, Clifford u2);

    template <typename T>
    Clifford clifford_inv(T op);

    std::tuple<Axis, Clifford> clifford_tconj(Clifford u);

    Pair<int> conj2(int a, int b);

    Tup4<int> conj3(int a, int b, int c);

    Tup4<int> cinv(int a, int b, int c);

    std::tuple<Axis, int, int> tconj(int a, int b);
}

#include "clifford.cpp"