#pragma once
#include <gmpxx.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <tuple>

using Integer = mpz_class;
using Rational = mpq_class;

namespace bmp = boost::multiprecision;

template <unsigned int N>
using Decimal = bmp::number<bmp::cpp_dec_float<N>, bmp::et_off>;

const int REAL_DIGITS = 100;

using Real = Decimal<REAL_DIGITS>;

template <typename T>
using Tup2 = std::tuple<T, T>;

template <typename T>
using Tup3 = std::tuple<T, T, T>;

template <typename T>
using Tup4 = std::tuple<T, T, T, T>;

template <typename T>
using Tup5 = std::tuple<T, T, T, T, T>;
