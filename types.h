#pragma once
#include <gmpxx.h>
#include <boost/multiprecision/cpp_dec_float.hpp>

using Integer = mpz_class; 
using Rational = mpq_class;

template <unsigned int N>
using Decimal = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<N>>;

const int REAL_DIGITS = 10;

using Real = Decimal<REAL_DIGITS>;
