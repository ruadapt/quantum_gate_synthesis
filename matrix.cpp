#include "matrix.h"
#include "ring.h"

namespace matrix
{
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

    template <typename T, size_t M, size_t N>
    Matrix<T, M, N> matrix_of_function(std::function<T(size_t, size_t)> f)
    {
        Matrix<T, M, N> result;
        for (size_t i = 0; i < M; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                result(i, j) = f(i, j);
            }
        }
        return result;
    }

    template <typename T, int N>
    Matrix<T, N, N> fromInteger(int x)
    {
        auto f = [=](size_t i, size_t j) -> T
        {
            if (i == j)
            {
                return ring::fromInteger<T>(x);
            }
            return 0;
        };
        return matrix_of_function<T, N, N>(f);
    }
}