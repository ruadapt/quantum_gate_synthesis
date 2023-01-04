#pragma once
#include "types.h"

template <typename T>
U2<T> matrix2x2(T x0, T x1, T x2, T x3)
{
    U2<T> m;
    m(0, 0) = x0;
    m(0, 1) = x1;
    m(1, 0) = x2;
    m(1, 1) = x3;
    return m;
}

#include "matrix.cpp"