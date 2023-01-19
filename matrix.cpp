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

    // TODO move to ring namespace
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

    // TODO find a more efficient way to do this with some kind of range/view
    template <typename T, int M, int N>
    List<T> get_row(Matrix<T, M, N> m, size_t r)
    {
        assert(r > 0 && r < M);
        List<T> row;
        row.reserve(N);
        for (size_t c = 0; c < N; c++)
        {
            row.push_back(m(r, c));
        }
        return row;
    }

    // TODO find a more efficient way to do this with some kind of range/view
    template <typename T, int M, int N>
    List<T> get_col(Matrix<T, M, N> m, size_t c)
    {
        assert(c > 0 && c < N);
        List<T> col;
        col.reserve(M);
        for (size_t r = 0; r < M; r++)
        {
            col.push_back(m(r, c));
        }
        return col;
    }

    /**
     * Break a matrix into its first column and the rest.
     */
    // TODO find a more efficient way to do this with some kind of range/view
    template <typename T, int M, int N>
    std::tuple<Matrix<T, M, 1>, Matrix<T, M, N - 1>> col_split(Matrix<T, M, N> m)
    {
        Matrix<T, M, 1> first_col;
        Matrix<T, M, N - 1> rest;
        for (size_t i = 0; i < M; i++)
        {
            first_col(i, 0) = m(i, 0);
        }
        for (size_t i = 0; i < M; i++)
        {
            for (size_t j = 1; j < N; j++)
            {
                rest(i, j - 1) = m(i, j);
            }
        }
        return {first_col, rest};
    }
}