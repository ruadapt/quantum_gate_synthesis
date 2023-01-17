#include "multiQubitSynthesis.h"
#include "matrix.h"
#include "ring.h"

namespace multi_qubit_synthesis
{
    namespace mat = matrix;

    Z2 residue(Integer n)
    {
        return ring::parity(n);
    }

    template <typename A, typename B>
    Omega<B> residue(Omega<A> o)
    {
        return Omega<B>(residue(o.A()), reside(o.b()), reside(o.c()), reside(o.d()));
    }

    template <typename A, typename B>
    Omega<B> residue(RootTwo<A> r)
    {
        return RootTwo<B>(residue(r.A()), reside(r.b()));
    }

    TwoLevel invert_twolevel(TwoLevel tl)
    {
        switch (tl.type())
        {
        case TL_X:
        {
            return make_TL_X(tl.i1(), tl.i2());
        }
        case TL_H:
        {
            return make_TL_H(tl.i1(), tl.i2());
        }
        case TL_T:
        {
            return make_TL_T(-tl.pow(), tl.i1(), tl.i2());
        }
        case TL_omega:
        {
            return make_TL_omega(-tl.pow(), tl.i1());
        }
        }
    }

    List<TwoLevel> invert_twolevels(List<TwoLevel> tls)
    {
        List<TwoLevel> inverted;
        std::transform(tls.begin(), tls.end(), std::back_inserter(inverted), invert_twolevel);
        std::reverse(inverted.begin(), inverted.end());
        return inverted;
    }

    template <typename T, int N>
    Matrix<T, N, N> twolevel_matrix(Pair<T> p1, Pair<T> p2, Index i, Index j)
    {
        T a, b, c, d;
        std::tie(a, b) = p1;
        std::tie(c, d) = p2;
        auto f = [=](size_t x, size_t y) -> T
        {
            if (x == i && y == i)
            {
                return a;
            }
            if (x == i && y == j)
            {
                return b;
            }
            if (x == j && y == i)
            {
                return c;
            }
            if (x == j && y == j)
            {
                return d;
            }
            if (x == y)
            {
                return T(1);
            }
            return T(0);
        };
        return mat::matrix_of_function<T, N, N>(f);
    }

    template <typename T, int N>
    Matrix<T, N, N> onelevel_matrix(T a, Index i)
    {
        auto f = [=](size_t x, size_t y) -> T
        {
            if (x == i && y == i)
            {
                return a;
            }
            if (x == y)
            {
                return T(1);
            }
            return T(0);
        };
        return mat::matrix_of_function<T, N, N>(f);
    }

    template <typename T, int N>
    Matrix<T, N, N> matrix_of_twolevel(TwoLevel tl)
    {
        switch (tl.type())
        {
        case TL_X:
        {
            return twolevel_matrix<T, N>(Pair<T>{0, 1}, Pair<T>{1, 0}, tl.i1(), tl.i2());
        }
        case TL_H:
        {
            T s = ring::rootHalf<T>();
            return twolevel_matrix<T, N>(Pair<T>{s, s}, Pair<T>{s, -s}, tl.i1(), tl.i2());
        }
        case TL_T:
        {
            T o = ring::powNonNeg(ring::omega<T>(), utils::mod(tl.pow(), 8));
            return twolevel_matrix<T, N>(Pair<T>{1, 0}, Pair<T>{0, o}, tl.i1(), tl.i2());
        }
        case TL_omega:
        {
            T o = ring::powNonNeg(ring::omega<T>(), utils::mod(tl.pow(), 8));
            return onelevel_matrix<T, N>(o, tl.i1());
        }
        }
    }

    template <typename T, int N>
    Matrix<T, N, N> matrix_of_twolevels(List<TwoLevel> gs)
    {
        Matrix<T, N, N> result = mat::fromInteger<T, N>(1);
        for (TwoLevel g : gs)
        {
            result = prod(result, matrix_of_twolevel<T, N>(g));
        }
        return result;
    }
}