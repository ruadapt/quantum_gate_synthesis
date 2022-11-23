#include "types.h"
#include "ring.h"
#include <iostream>
#include <gmpxx.h>
#include <cmath>
#include <boost/multiprecision/number.hpp>

namespace bp = boost::multiprecision;

template <int N>
std::string toString(Decimal<N> d)
{
    return d.str(N);
}

int main()
{
    Decimal<50> x("1.2345678901122334455667788990012341234123412341234123412341234");
    std::cout << "digits: " << std::numeric_limits<Decimal<50>>::digits10 << std::endl;
    // std::cout << std::setprecision(std::numeric_limits<Decimal<50>>::digits10) << x << std::endl;
    std::cout << "x = " << toString<50>(x) << std::endl;

    Real y = 1.5;
    auto s1 = bp::pow(x, 2);
    Real s2 = pow(y, 2);
    std::cout << "x^2 = " << toString<50>(s1) << std::endl;
    std::cout << "y^2 = " << s2 << std::endl;
}