#pragma once
#include "types.h"

class Clifford
{
public:
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
        return "Clifford(" + std::to_string(a()) + ", " + std::to_string(b()) + ", " + std::to_string(c()) + ", " + std::to_string(d()) + ")";
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

#include "clifford.cpp"