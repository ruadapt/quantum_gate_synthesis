#pragma once
#include <gmpxx.h>
#include <boost/multiprecision/cpp_dec_float.hpp>

using Integer = mpz_class;
using Rational = mpq_class;

namespace bmp = boost::multiprecision;

template <unsigned int N>
using Decimal = bmp::number<bmp::cpp_dec_float<N>, bmp::et_off>;

const int REAL_DIGITS = 100;

using Real = Decimal<REAL_DIGITS>;
