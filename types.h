#pragma once
#include <gmpxx.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <tuple>
#include <variant>
#include <optional>

using Integer = mpz_class;
using Rational = mpq_class;

namespace bmp = boost::multiprecision;

template <unsigned int N>
using Decimal = bmp::number<bmp::cpp_dec_float<N>, bmp::et_off>;

const int REAL_DIGITS = 100;

using Real = Decimal<REAL_DIGITS>;

template <typename T>
using List = std::vector<T>;

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

template <typename A, typename B>
using Either = std::variant<std::monostate, A, B>;

template <typename A, typename B>
A fst(Either<A, B> e) { return std::get<1>(e); }

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
using Tuple2By2 = std::tuple<std::tuple<T, T>, std::tuple<T, T>>;

template <typename T>
using Operator = Matrix<T, 2, 2>;

template <typename T>
using OperatorPair = std::tuple<Operator<T>, Operator<T>>;