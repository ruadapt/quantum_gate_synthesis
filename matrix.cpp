#include "matrix.h"

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