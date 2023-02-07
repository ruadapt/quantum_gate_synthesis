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
        return Omega<B>(residue(o.a()), residue(o.b()), residue(o.c()), residue(o.d()));
    }

    template <typename A, typename B>
    RootTwo<B> residue(RootTwo<A> r)
    {
        return RootTwo<B>(residue(r.A()), residue(r.b()));
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

    template <typename T, size_t N>
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

    template <typename T, size_t N>
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

    template <typename T, size_t N>
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
            T s = ring::roothalf<T>();
            return twolevel_matrix<T, N>(Pair<T>{s, s}, Pair<T>{s, -s}, tl.i1(), tl.i2());
        }
        case TL_T:
        {
            T o = ring::pow_non_neg(ring::omega<T>(), utils::mod(tl.pow(), 8));
            return twolevel_matrix<T, N>(Pair<T>{1, 0}, Pair<T>{0, o}, tl.i1(), tl.i2());
        }
        case TL_omega:
        {
            T o = ring::pow_non_neg(ring::omega<T>(), utils::mod(tl.pow(), 8));
            return onelevel_matrix<T, N>(o, tl.i1());
        }
        }
    }

    template <typename T, size_t N>
    Matrix<T, N, N> matrix_of_twolevels(List<TwoLevel> gs)
    {
        Matrix<T, N, N> result = mat::fromInteger<T, N>(1);
        for (TwoLevel g : gs)
        {
            result = prod(result, matrix_of_twolevel<T, N>(g));
        }
        return result;
    }

    template <typename T>
    List<T> list_insert(Index n, T x, List<T> lst)
    {
        if (n >= lst.size())
        {
            return lst;
        }
        List<T> result = lst;
        result.at(n) = x;
        return result;
    }

    template <typename T>
    List<T> transform_at(std::function<T(T)> op, Index i, List<T> lst)
    {
        if (i >= lst.size())
        {
            return lst;
        }
        T x = lst.at(i);
        List<T> result = lst;
        result.at(i) = op(x);
        return result;
    }

    template <typename T>
    List<T> transform_at2(std::function<Pair<T>(Pair<T>)> op, Index i, Index j, List<T> lst)
    {
        if (i >= lst.size())
        {
            return lst;
        }
        T x = lst.at(i);
        T y = lst.at(j);
        T x_prime, y_prime;
        std::tie(x_prime, y_prime) = op(Pair<T>{x, y});
        List<T> result = lst;
        result.at(i) = x_prime;
        result.at(j) = y_prime;
        return result;
    }

    template <typename T>
    std::tuple<List<Pair<T>>, Maybe<T>> list_pairs(List<T> lst)
    {
        List<Pair<T>> pairs;
        bool even = lst.size() % 2 == 0;
        size_t n = even ? lst.size() : (lst.size() - 1);
        for (size_t i = 0; i < n; i += 2)
        {
            pairs.push_back(Pair<T>{lst.at(i), lst.at(i + 1)});
        }
        Maybe<T> last = even ? Maybe<T>() : Maybe<T>(lst.at(lst.size() - 1));
        return {pairs, last};
    }

    Maybe<int> log_omega(ZOmega z)
    {
        if (z == ZOmega(0, 0, 0, 1))
        {
            return 0;
        }
        if (z == ZOmega(0, 0, 1, 0))
        {
            return 1;
        }
        if (z == ZOmega(0, 1, 0, 0))
        {
            return 2;
        }
        if (z == ZOmega(1, 0, 0, 0))
        {
            return 3;
        }
        if (z == ZOmega(0, 0, 0, -1))
        {
            return 4;
        }
        if (z == ZOmega(0, 0, -1, 0))
        {
            return 5;
        }
        if (z == ZOmega(0, -1, 0, 0))
        {
            return 6;
        }
        if (z == ZOmega(-1, 0, 0, 0))
        {
            return 7;
        }
        return Maybe<int>();
    }

    template <typename T>
    T omega_power(int n, T x)
    {
        return x * ring::pow_non_neg(ring::omega<T>(), utils::mod(n, 8));
    }

    ZOmega reduce_ZOmega(ZOmega z)
    {
        Integer a, b, c, d;
        a = z.a();
        b = z.b();
        c = z.c();
        d = z.d();
        if (ring::even(a - c) && ring::even(b - d))
        {
            Integer a_prime = utils::div(b - d, 2_mpz);
            Integer b_prime = utils::div(c + a, 2_mpz);
            Integer c_prime = utils::div(b + d, 2_mpz);
            Integer d_prime = utils::div(c - a, 2_mpz);
            return ZOmega(a_prime, b_prime, c_prime, d_prime);
        }
        else
        {
            throw std::invalid_argument("Argument is not reducible");
        }
    }

    Pair<ZOmega> opX_zomega(Pair<ZOmega> p)
    {
        return Pair<ZOmega>(snd(p), fst(p));
    }

    Pair<ZOmega> opH_zomega(Pair<ZOmega> p)
    {
        ZOmega x, y;
        std::tie(x, y) = p;
        return Pair<ZOmega>{reduce_ZOmega(x + y), reduce_ZOmega(x - y)};
    }

    List<ZOmega> apply_twolevel_zomega(TwoLevel tl, List<ZOmega> w)
    {
        switch (tl.type())
        {
        case TL_X:
        {
            return transform_at2<ZOmega>(opX_zomega, tl.i1(), tl.i2(), w);
        }
        case TL_H:
        {
            return transform_at2<ZOmega>(opH_zomega, tl.i1(), tl.i2(), w);
        }
        case TL_T:
        {
            int k = tl.pow();
            auto f = [=](ZOmega z) -> ZOmega
            {
                return omega_power(k, z);
            };
            return transform_at<ZOmega>(f, tl.i2(), w);
        }
        case TL_omega:
        {
            int k = tl.pow();
            auto f = [=](ZOmega z) -> ZOmega
            {
                return omega_power(k, z);
            };
            return transform_at<ZOmega>(f, tl.i1(), w);
        }
        }
    }

    List<ZOmega> apply_twolevels_zomega(List<TwoLevel> gs, List<ZOmega> w)
    {
        return utils::foldr<TwoLevel, List<ZOmega>>(apply_twolevel_zomega, w, gs);
    }

    std::tuple<ResidueType, int> residue_type_shift(Omega<Z2> r)
    {
        if (r == Omega<Z2>(0, 0, 0, 0))
        {
            return {RT_0000, 0};
        }
        if (r == Omega<Z2>(0, 0, 0, 1))
        {
            return {RT_0001, 0};
        }
        if (r == Omega<Z2>(0, 0, 1, 0))
        {
            return {RT_0001, 1};
        }
        if (r == Omega<Z2>(0, 0, 1, 1))
        {
            return {RT_1010, 0};
        }
        if (r == Omega<Z2>(0, 1, 0, 0))
        {
            return {RT_0001, 2};
        }
        if (r == Omega<Z2>(0, 1, 0, 1))
        {
            return {RT_0000, 0};
        }
        if (r == Omega<Z2>(0, 1, 1, 0))
        {
            return {RT_1010, 1};
        }
        if (r == Omega<Z2>(0, 1, 1, 1))
        {
            return {RT_0001, 3};
        }
        if (r == Omega<Z2>(1, 0, 0, 0))
        {
            return {RT_0001, 3};
        }
        if (r == Omega<Z2>(1, 0, 0, 1))
        {
            return {RT_1010, 3};
        }
        if (r == Omega<Z2>(1, 0, 1, 0))
        {
            return {RT_0000, 0};
        }
        if (r == Omega<Z2>(1, 0, 1, 1))
        {
            return {RT_0001, 2};
        }
        if (r == Omega<Z2>(1, 1, 0, 0))
        {
            return {RT_1010, 2};
        }
        if (r == Omega<Z2>(1, 1, 0, 1))
        {
            return {RT_0001, 1};
        }
        if (r == Omega<Z2>(1, 1, 1, 0))
        {
            return {RT_0001, 0};
        }
        if (r == Omega<Z2>(1, 1, 1, 1))
        {
            return {RT_0000, 0};
        }
        // Since Z2 can only be 0 or 1, we should never reach this point.
        throw std::invalid_argument("This exception shouldn't be reachable");
    }

    ResidueType residue_type(Omega<Z2> r)
    {
        return fst(residue_type_shift(r));
    }

    int residue_shift(Omega<Z2> r)
    {
        return snd(residue_type_shift(r));
    }

    int residue_offset(Omega<Z2> a, Omega<Z2> b)
    {
        return utils::mod(residue_shift(a) - residue_shift(b), 4);
    }

    bool reducible(Omega<Z2> r)
    {
        return (r.a() == r.c()) && (r.b() == r.d());
    }

    List<TwoLevel> row_step(Pair<std::tuple<Index, Omega<Z2>, ZOmega>> p)
    {
        std::tuple<Index, Omega<Z2>, ZOmega> t1, t2;
        std::tie(t1, t2) = p;
        Index i, j;
        Omega<Z2> a, b;
        ZOmega x, y;
        std::tie(i, a, x) = t1;
        std::tie(j, b, y) = t2;
        if (reducible(a) && reducible(b))
        {
            return List<TwoLevel>{};
        }
        int offs = residue_offset(b, a);
        if (offs != 0)
        {
            ZOmega y_prime = omega_power(-offs, y);
            Omega<Z2> b_prime = residue<Integer, Z2>(y_prime);
            std::tuple<Index, Omega<Z2>, ZOmega> t2_prime{j, b_prime, y_prime};
            return utils::cons(make_TL_T(offs, i, j), row_step({t1, t2_prime}));
        }
        ZOmega x1, y1;
        std::tie(x1, y1) = opH_zomega(Pair<ZOmega>{x, y});
        Omega<Z2> a1 = residue<Integer, Z2>(x1);
        Omega<Z2> b1 = residue<Integer, Z2>(y1);
        std::tuple<Index, Omega<Z2>, ZOmega> t1_prime{i, a1, x1};
        std::tuple<Index, Omega<Z2>, ZOmega> t2_prime{j, b1, y1};
        return utils::cons(make_TL_H(i, j), row_step({t1_prime, t2_prime}));
    }

    List<TwoLevel> reduce_column_aux(List<ZOmega> w, Integer k, Index i)
    {
        if (k == 0)
        {
            int non_zero_count = 0;
            Index j;
            for (Index idx = 0; idx < w.size(); idx++)
            {
                if (w.at(idx) != 0)
                {
                    j = idx;
                    non_zero_count++;
                    if (non_zero_count > 1)
                    {
                        break;
                    }
                }
            }
            if (non_zero_count != 1)
            {
                throw std::invalid_argument("Input was not a unit vector");
            }
            ZOmega wj = w.at(j);
            Maybe<int> log = log_omega(wj);
            if (!log.has_value())
            {
                throw std::invalid_argument("Input was not a unit vector");
            }
            int l = log.value();
            List<TwoLevel> m1 = (i == j) ? List<TwoLevel>{} : List<TwoLevel>{make_TL_X(i, j)};
            List<TwoLevel> m2{make_TL_omega(l, i)};
            return utils::concat(m1, m2);
        }
        auto residue_omega = [](ZOmega z) -> Omega<Z2>
        { return residue<Integer, Z2>(z); };
        List<Omega<Z2>> res = utils::map<ZOmega, Omega<Z2>>(residue_omega, w);
        // Make a temporary type alias since this type is used several times.
        using C = std::tuple<Index, Omega<Z2>, ZOmega>;
        List<C> res1010;
        List<C> res0001;
        for (Index i = 0; i < res.size(); i++)
        {
            Omega<Z2> a = res.at(i);
            ZOmega x = w.at(i);
            ResidueType rt = residue_type(a);
            if (rt == RT_1010)
            {
                res1010.push_back({i, a, x});
            }
            else if (rt == RT_0001)
            {
                res0001.push_back({i, a, x});
            }
        }
        List<Pair<C>> res1010_pairs, res0001_pairs;
        Maybe<C> extra1010, extra0001;
        std::tie(res1010_pairs, extra1010) = list_pairs(res1010);
        std::tie(res0001_pairs, extra0001) = list_pairs(res0001);
        if (extra1010.has_value() || extra0001.has_value())
        {
            throw std::invalid_argument("Input was not a unit vector");
        }
        List<TwoLevel> m1010 = utils::concat(utils::map<Pair<C>, List<TwoLevel>>(row_step, res1010_pairs));
        List<TwoLevel> m0001 = utils::concat(utils::map<Pair<C>, List<TwoLevel>>(row_step, res0001_pairs));
        List<TwoLevel> gates = utils::concat(m1010, m0001);
        List<ZOmega> applied = apply_twolevels_zomega(invert_twolevels(gates), w);
        List<ZOmega> w_prime = utils::map<ZOmega, ZOmega>(reduce_ZOmega, applied);
        return utils::concat(gates, reduce_column_aux(w_prime, k - 1, i));
    }

    std::tuple<List<ZOmega>, Integer> denomexp_decompose(List<DOmega> lst)
    {
        List<Integer> exponents = utils::map<DOmega, Integer>(ring::denomexp<DOmega>, lst);
        Integer k = utils::max(exponents, 0_mpz);
        auto factor = [k](DOmega z) -> DOmega
        { return ring::denomexp_factor<DOmega>(z, k); };
        List<DOmega> factored = utils::map<DOmega, DOmega>(factor, lst);
        List<ZOmega> whole = utils::map<DOmega, ZOmega>(ring::to_whole<DOmega, ZOmega>, factored);
        return {whole, k};
    }

    template <size_t N>
    List<TwoLevel> reduce_column(Matrix<DOmega, N, 1> v, Index i)
    {
        List<DOmega> vlist = mat::get_col(v, 0);
        List<ZOmega> w;
        Integer k;
        std::tie(w, k) = denomexp_decompose(vlist);
        return reduce_column_aux(w, k, i);
    }

    template <size_t M, size_t N>
    List<TwoLevel> synthesis_nqubit_aux(Matrix<DOmega, M, N> m, Index i)
    {
        // We need the constexpr because otherwise, N - 1 will cause an overflow
        // since size_t is unsigned, which causes an error when allocating cs.
        if constexpr (N == 0)
        {
            return List<TwoLevel>{};
        }
        else
        {
            Matrix<DOmega, M, 1> c;
            Matrix<DOmega, M, N - 1> cs;
            std::tie(c, cs) = mat::col_split(m);
            List<TwoLevel> gates = reduce_column(c, i);
            Matrix<DOmega, M, M> gates_matrix = matrix_of_twolevels<DOmega, M>(invert_twolevels(gates));
            Matrix<DOmega, M, N - 1> m_prime = prod(gates_matrix, cs);
            return utils::concat(gates, synthesis_nqubit_aux(m_prime, i + 1));
        }
    }

    template <size_t N>
    List<TwoLevel> synthesis_nqubit(Matrix<DOmega, N, N> m)
    {
        return synthesis_nqubit_aux(m, 0);
    }
}