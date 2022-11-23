#pragma once
#include <gmpxx.h>
#include <boost/multiprecision/cpp_dec_float.hpp>

using Integer = mpz_class; 
using Rational = mpq_class;

template <int digits>
using Decimal = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<digits>>;

using Real = double;