/** \file types.h
 */
#pragma once
#include <gmpxx.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <tuple>
#include <variant>
#include <optional>

#ifndef REAL_DIGITS
    #define REAL_DIGITS 100
#endif

/**
 * @brief Arbitrary precision integer type.
 */
using Integer = mpz_class;

/**
 * @brief Arbitrary precision rational type.
 */
using Rational = mpq_class;

namespace bmp = boost::multiprecision;

/**
 * @brief A real number with N digits of precision.
 */
template <unsigned int N>
using Decimal = bmp::number<bmp::cpp_dec_float<N>, bmp::et_off>;

/**
 * The default real number type, with REAL_DIGITS digits of precision.
 */
using Real = Decimal<REAL_DIGITS>;

template <typename T>
using List = std::vector<T>;

template <typename T>
struct is_list : std::false_type
{
};

template <typename T>
struct is_list<List<T>> : std::true_type
{
};

template <typename T> constexpr bool is_list_v = is_list<T>::value;

template <typename T>
std::ostream &operator<<(std::ostream &os, const List<T> &l)
{
    os << "[";
    for (size_t i = 0; i < l.size(); i++)
    {
        os << l.at(i);
        if (i != l.size() - 1)
        {
            os << ", ";
        }
    }
    os << "]";
    return os;
}

template <typename T>
using Maybe = std::optional<T>;

template <typename T>
std::ostream &operator<<(std::ostream &os, const Maybe<T> &m)
{
    if (m.has_value())
    {
        os << "Just(" << m.value() << ")";
    }
    else
    {
        os << "Nothing";
    }
    return os;
}

template <typename T>
using Pair = std::tuple<T, T>;

template <typename A, typename B>
A fst(std::tuple<A, B> t) { return std::get<0>(t); }

template <typename A, typename B>
B snd(std::tuple<A, B> t) { return std::get<1>(t); }

template <typename T>
std::ostream &operator<<(std::ostream &os, const Pair<T> &p)
{
    os << "(" << fst(p) << ", " << snd(p) << ")";
    return os;
}

template <typename T>
using Tup3 = std::tuple<T, T, T>;

template <typename T>
using Tup4 = std::tuple<T, T, T, T>;

template <typename T>
using Tup5 = std::tuple<T, T, T, T, T>;

/**
 * std::monostate is there so that this type is still default-constructible even if
 * A and B aren't.
 */
template <typename A, typename B>
using Either = std::variant<std::monostate, A, B>;

/**
 * @brief Get the A value.
 */
template <typename A, typename B>
A fst(Either<A, B> e) { return std::get<1>(e); }

/**
 * @brief Get the B value.
 */
template <typename A, typename B>
B snd(Either<A, B> e) { return std::get<2>(e); }

template <typename T, size_t N>
using Vector = boost::numeric::ublas::c_vector<T, N>;

template <typename T, size_t M, size_t N>
using Matrix = boost::numeric::ublas::c_matrix<T, M, N>;

template <typename T>
using U2 = Matrix<T, 2, 2>;

template <typename T>
using SO3 = Matrix<T, 3, 3>;

template <typename T>
using Point = std::tuple<T, T>;

template <typename T>
using Tuple2By2 = std::tuple<Pair<T>, Pair<T>>;

template <typename T>
using Operator = Matrix<T, 2, 2>;

template <typename T>
using OperatorPair = std::tuple<Operator<T>, Operator<T>>;