#pragma once
#include <boost/numeric/ublas/matrix.hpp>

template <typename T, int N>
class Vector: public boost::numeric::ublas::c_vector<T, N>
{
};

template <typename T, int M, int N>
class Matrix: public boost::numeric::ublas::c_matrix<T, M, N>
{
};