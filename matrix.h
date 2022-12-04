#pragma once
#include <boost/numeric/ublas/matrix.hpp>

template <typename T, int N>
using Vector = boost::numeric::ublas::c_vector<T, N>;

template <typename T, int M, int N>
using Matrix = boost::numeric::ublas::c_matrix<T, M, N>;