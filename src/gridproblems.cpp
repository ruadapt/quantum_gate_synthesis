#include "ring.h"
#include "quadratic.h"
#include <boost/math/constants/constants.hpp>
#include <vector>

namespace constants = boost::math::constants;
namespace bmp = boost::multiprecision;

namespace gridprob
{
    namespace mat = matrix;

    template <typename T>
    Operator<T> operator_pow(Operator<T> base, int exp)
    {
        assert(exp >= 0);
        Operator<T> result;
        if (exp == 0)
        {
            return mat::matrix2x2<T>(1, 0, 0, 1);
        }
        if (exp == 1)
        {
            return base;
        }
        Operator<T> new_base = prod(base, base);
        if (exp % 2 == 0)
        {
            int new_exp = exp / 2;
            return operator_pow(new_base, new_exp);
        }
        int new_exp = (exp - 1) / 2;
        return prod(base, operator_pow(new_base, new_exp));
    }

    // TODO make it possible to use power methods in ring on operators
    template <typename T>
    Operator<T> operator_pow(Operator<T> base, Integer exp)
    {
        return operator_pow(base, utils::to_int(exp));
    }

    // TODO move to ring namespace
    /**
     * Note that this isn't just an elementwise adjoint. There's also a transpose operation afterwards,
     * which is why the assignments to result(1, 0) and result(0, 1) are flipped. This behavior
     * matches what's done in Haskell.
     */
    template <typename T>
    Operator<T> adj(Operator<T> op)
    {
        Operator<T> result;
        result(0, 0) = ring::adj(op(0, 0));
        result(1, 0) = ring::adj(op(0, 1));
        result(0, 1) = ring::adj(op(1, 0));
        result(1, 1) = ring::adj(op(1, 1));
        return result;
    }

    // TODO move to ring namespace
    template <typename T>
    Operator<T> adj2(Operator<T> op)
    {
        Operator<T> result;
        result(0, 0) = ring::adj2<T>(op(0, 0));
        result(0, 1) = ring::adj2<T>(op(0, 1));
        result(1, 0) = ring::adj2<T>(op(1, 0));
        result(1, 1) = ring::adj2<T>(op(1, 1));
        return result;
    }

    template <typename T>
    List<T> adj(List<T> v)
    {
        List<T> result(v.size());
        std::transform(v.begin(), v.end(), result.begin(), [](T x)
                       { return ring::adj<T>(x); });
        return result;
    }

    template <typename T>
    List<T> adj2(List<T> v)
    {
        List<T> result(v.size());
        std::transform(v.begin(), v.end(), result.begin(), [](T x)
                       { return ring::adj2<T>(x); });
        return result;
    }

    double logBase(double b, double x)
    {
        return log(x) / log(b);
    }

    template <typename T>
    T lambda()
    {
        return T(1) + ring::roottwo<T>();
    }

    template <typename T>
    T lambda_inv()
    {
        return ring::roottwo<T>() - 1;
    }

    template <typename T>
    T lambdapower(Integer k)
    {
        if (k >= 0)
        {
            return ring::pow_non_neg(lambda<T>(), k);
        }
        return ring::pow_non_neg(lambda_inv<T>(), -k);
    }

    template <typename T>
    T signpower(Integer k)
    {
        return ring::even(k) ? T(1) : T(-1);
    }

    template <typename T>
    bool within(const T x, const T low, const T high)
    {
        return (low <= x && x <= high);
    }

    template <typename T>
    Pair<T> fatten_interval(Pair<T> interval)
    {
        T x = std::get<0>(interval);
        T y = std::get<1>(interval);
        T epsilon = 0.0001 * (y - x);
        return std::make_tuple(x - epsilon, y + epsilon);
    }

    template <typename T>
    std::tuple<Integer, T> floorlog(T b, T x)
    {
        if (x <= T(0))
        {
            throw std::invalid_argument("x must be greater than 0");
        }
        if (T(1) <= x && x < b)
        {
            return std::make_tuple(0, x);
        }
        if (T(1) <= x * b && x < T(1))
        {
            return std::make_tuple(-1, b * x);
        }
        Integer n;
        T r;
        std::tie(n, r) = floorlog<T>(b * b, x);
        if (r < b)
        {
            return std::make_tuple(2 * n, r);
        }
        return std::make_tuple(2 * n + 1, r * ring::recip<T>(b));
    }

    template <typename T>
    double logBase_double(T b, T x)
    {
        if (b > 1)
        {
            Integer n;
            T r;
            std::tie(n, r) = floorlog(b, x);
            // TODO check the conversions here
            return ring::fromInteger<double>(n) + logBase(double(b), double(r));
        }
        if (b <= 0)
        {
            return NAN;
        }
        if (b == 1)
        {
            return std::numeric_limits<double>::infinity();
        }
        return -logBase_double(ring::recip(b), x);
    }

    template <typename T>
    T iprod(Point<T> p1, Point<T> p2)
    {
        return fst(p1) * fst(p2) + snd(p1) * snd(p2);
    }

    template <typename T>
    Operator<T> special_inverse(Operator<T> opG)
    {
        T a, b, c, d;
        a = opG(0, 0);
        b = opG(0, 1);
        c = opG(1, 0);
        d = opG(1, 1);
        return det(opG) * mat::matrix2x2(d, -b, -c, a);
    }

    template <typename T>
    List<ZRootTwo> gridpoints_internal(T x0, T x1, T y0, T y1)
    {
        T dx = x1 - x0;
        T dy = y1 - y0;

        if (dy < 0 || dx < 0)
        {
            return List<ZRootTwo>{};
        }

        auto base_case = [=]
        {
            Integer amin = ring::ceiling_of((x0 + y0) * ring::recip<T>(2));
            Integer amax = ring::floor_of((x1 + y1) * ring::recip<T>(2));
            List<ZRootTwo> results;
            for (Integer a = amin; a <= amax; a++)
            {
                T a_T = ring::fromInteger<T>(a);
                Integer bmin = ring::ceiling_of((a_T - y1) * ring::recip(ring::roottwo<T>()));
                Integer bmax = ring::floor_of((a_T - y0) * ring::recip(ring::roottwo<T>()));
                for (Integer b = bmin; b <= bmax; b++)
                {
                    results.push_back(ZRootTwo(a, b));
                }
            }
            return results;
        };

        if (dy == 0 && dx == 0)
        {
            return base_case();
        }

        if (dy == 0 && dx > 0)
        {
            List<ZRootTwo> sub_results = gridpoints_internal(y0, y1, x0, x1);
            List<ZRootTwo> results(sub_results.size());
            std::transform(sub_results.begin(), sub_results.end(), results.begin(), [](ZRootTwo z)
                           { return z.adj2(); });
            return results;
        }

        // floorlog requires a positive argument, so we don't compute it until
        // we've ruled out the possibility that dy <= 0.
        int n = utils::to_int(std::get<0>(floorlog<T>(lambda<T>(), dy)));
        int m = -n;

        T lambda_T = lambda<T>();
        T lambda_inv_T = lambda_inv<T>();
        T lambda_n_T = ring::pow_int(lambda_T, n);
        T lambda_bul_n_T = ring::pow_int(-lambda_inv_T, n);
        T lambda_inv_m_T = ring::pow_int(lambda_inv_T, m);
        T lambda_bul_inv_m_T = ring::pow_int(-lambda_T, m);
        ZRootTwo lambda_Z = lambda<ZRootTwo>();
        ZRootTwo lambda_inv_Z = lambda_inv<ZRootTwo>();

        if (dy >= lambda_T && ring::even(n))
        {
            List<ZRootTwo> sub_results = gridpoints_internal(
                lambda_n_T * x0, lambda_n_T * x1, lambda_bul_n_T * y0, lambda_bul_n_T * y1);
            List<ZRootTwo> results(sub_results.size());
            ZRootTwo lambda_inv_n_Z = ring::pow_non_neg(lambda_inv_Z, n);
            std::transform(sub_results.begin(), sub_results.end(), results.begin(),
                           [=](ZRootTwo z)
                           { return lambda_inv_n_Z * z; });
            return results;
        }

        if (dy >= lambda_T && ring::odd(n))
        {
            List<ZRootTwo> sub_results = gridpoints_internal(
                lambda_n_T * x0, lambda_n_T * x1, lambda_bul_n_T * y1, lambda_bul_n_T * y0);
            List<ZRootTwo> results(sub_results.size());
            ZRootTwo lambda_inv_n_Z = ring::pow_non_neg(lambda_inv_Z, n);
            std::transform(sub_results.begin(), sub_results.end(), results.begin(),
                           [=](ZRootTwo z)
                           { return lambda_inv_n_Z * z; });
            return results;
        }

        if (dy > 0 && dy < 1 && ring::even(n))
        {
            List<ZRootTwo> sub_results = gridpoints_internal(
                lambda_inv_m_T * x0, lambda_inv_m_T * x1, lambda_bul_inv_m_T * y0, lambda_bul_inv_m_T * y1);
            List<ZRootTwo> results(sub_results.size());
            ZRootTwo lambda_m_Z = ring::pow_non_neg(lambda_Z, m);
            std::transform(sub_results.begin(), sub_results.end(), results.begin(),
                           [=](ZRootTwo z)
                           { return lambda_m_Z * z; });
            return results;
        }

        if (dy > 0 && dy < 1 && ring::odd(n))
        {
            List<ZRootTwo> sub_results = gridpoints_internal(
                lambda_inv_m_T * x0, lambda_inv_m_T * x1, lambda_bul_inv_m_T * y1, lambda_bul_inv_m_T * y0);
            List<ZRootTwo> results(sub_results.size());
            ZRootTwo lambda_m_Z = ring::pow_non_neg(lambda_Z, m);
            std::transform(sub_results.begin(), sub_results.end(), results.begin(),
                           [=](ZRootTwo z)
                           { return lambda_m_Z * z; });
            return results;
        }

        return base_case();
    }

    template <typename T>
    List<ZRootTwo> gridpoints(T x0, T x1, T y0, T y1)
    {
        Integer a = ring::div(ring::floor_of(x0 + y0), 2_mpz);
        Integer b = ring::div(ring::floor_of(ring::roottwo<T>() * (x0 - y0)), 4_mpz);
        ZRootTwo alpha = ZRootTwo(a, b);
        T xoff = ring::fromZRootTwo<T>(alpha);
        T yoff = ring::fromZRootTwo<T>(ring::adj2(alpha));
        T x0new = x0 - xoff;
        T x1new = x1 - xoff;
        T y0new = y0 - yoff;
        T y1new = y1 - yoff;

        List<ZRootTwo> raw_points = gridpoints_internal(x0new, x1new, y0new, y1new);
        List<ZRootTwo> shifted(raw_points.size());
        std::transform(raw_points.begin(), raw_points.end(), shifted.begin(), [=](ZRootTwo z)
                       { return z + alpha; });
        List<ZRootTwo> results;
        std::copy_if(shifted.begin(), shifted.end(), std::back_inserter(results), [=](ZRootTwo z)
                     {
                        T z1 = ring::fromZRootTwo<T>(z);
                        T z2 = ring::fromZRootTwo<T>(z.adj2());
                        return within<T>(z1, x0, x1) && within<T>(z2, y0, y1); });
        return results;
    }

    /**
     * @brief A convenient version of the other gridpoints_scaled that takes in two
     * pairs instead of 4 values.
     */
    template <typename T>
    List<DRootTwo> gridpoints_scaled(T x0, T x1, T y0, T y1, Integer k)
    {
        if (k < 0)
        {
            throw std::invalid_argument("gridpoints_scaled: k >= 0 is required");
        }
        T scale_inv_T = ring::pow_non_neg(ring::roottwo<T>(), k);
        DRootTwo scale_D = ring::pow_non_neg(ring::roothalf<DRootTwo>(), k);
        DRootTwo scale_inv_D = ring::pow_non_neg(ring::roottwo<DRootTwo>(), k);
        T x0new = scale_inv_T * x0;
        T x1new = scale_inv_T * x1;
        T y0new = ring::even(k) ? (scale_inv_T * y0) : (-scale_inv_T * y1);
        T y1new = ring::even(k) ? (scale_inv_T * y1) : (-scale_inv_T * y0);
        List<ZRootTwo> w = gridpoints(x0new, x1new, y0new, y1new);
        List<DRootTwo> results(w.size());
        std::transform(w.begin(), w.end(), results.begin(), [=](ZRootTwo z)
                       { return scale_D * ring::fromZRootTwo<DRootTwo>(z); });
        return results;
    }

    template <typename T>
    List<DRootTwo> gridpoints_scaled(Pair<T> intX, Pair<T> intY, Integer k)
    {
        T x0 = std::get<0>(intX);
        T x1 = std::get<1>(intX);
        T y0 = std::get<0>(intY);
        T y1 = std::get<1>(intY);
        return gridpoints_scaled(x0, x1, y0, y1, k);
    }

    template <typename T>
    List<DRootTwo> gridpoints_scaled_parity(DRootTwo beta, T x0, T x1, T y0, T y1, Integer k)
    {
        if (k < 1)
        {
            throw std::invalid_argument("gridpoints_scaled_parity: k >= 1 is required");
        }
        if (ring::denomexp(beta) <= k - 1)
        {
            return gridpoints_scaled(x0, x1, y0, y1, k - 1);
        }
        DRootTwo offs = ring::pow_non_neg(ring::roothalf<DRootTwo>(), k);
        DRootTwo offs_bul = ring::adj2(offs);
        T offs_prime = ring::fromDRootTwo<T>(offs);
        T offs_bul_prime = ring::fromDRootTwo<T>(offs_bul);
        List<DRootTwo> z = gridpoints_scaled(x0 + offs_prime, x1 + offs_prime, y0 + offs_bul_prime, y1 + offs_bul_prime, k - 1);
        List<DRootTwo> results(z.size());
        std::transform(z.begin(), z.end(), results.begin(), [=](DRootTwo d)
                       { return d - offs; });
        return results;
    }

    template <typename T>
    Point<T> point_fromDRootTwo(Point<DRootTwo> p)
    {
        return std::make_tuple(ring::fromDRootTwo<T>(fst(p)), ring::fromDRootTwo<T>(snd(p)));
    }

    template <typename T>
    ConvexSet<T> unitdisk()
    {
        Operator<T> op = mat::matrix2x2<T>(1, 0, 0, 1);
        Point<T> p = std::make_tuple<T, T>(0, 0);
        Ellipse<T> el = Ellipse<T>(op, p);

        std::function<Maybe<Pair<T>>(Point<DRootTwo>, Point<DRootTwo>)> intersect = [](Point<DRootTwo> p, Point<DRootTwo> v)
        {
            QRootTwo a = ring::fromDRootTwo<QRootTwo>(iprod<DRootTwo>(v, v));
            QRootTwo b = ring::fromDRootTwo<QRootTwo>(DRootTwo(2) * iprod<DRootTwo>(v, p));
            QRootTwo c = ring::fromDRootTwo<QRootTwo>(iprod<DRootTwo>(p, p) - 1);
            Maybe<Pair<T>> q = quadratic<QRootTwo>(a, b, c);
            return q;
        };

        std::function<bool(Point<DRootTwo>)> test = [](Point<DRootTwo> p)
        {
            DRootTwo x, y;
            std::tie(x, y) = p;
            return x * x + y * y <= 1;
        };

        return ConvexSet<T>(el, test, intersect);
    }

    template <typename T>
    ConvexSet<T> disk(DRootTwo s)
    {
        if (s <= 0)
        {
            throw std::invalid_argument("s > 0 is required");
        }
        T r = ring::recip<T>(ring::fromDRootTwo<T>(s));
        Operator<T> op = mat::matrix2x2<T>(r, 0, 0, r);
        Point<T> p = std::make_tuple<T, T>(0, 0);
        Ellipse<T> el = Ellipse<T>(op, p);

        std::function<Maybe<Pair<T>>(Point<DRootTwo>, Point<DRootTwo>)> intersect = [=](Point<DRootTwo> p, Point<DRootTwo> v)
        {
            QRootTwo a = ring::fromDRootTwo<QRootTwo>(iprod<DRootTwo>(v, v));
            QRootTwo b = ring::fromDRootTwo<QRootTwo>(DRootTwo(2) * iprod<DRootTwo>(v, p));
            QRootTwo c = ring::fromDRootTwo<QRootTwo>(iprod<DRootTwo>(p, p) - s);
            Maybe<std::tuple<Real, Real>> q = quadratic<QRootTwo>(a, b, c);
            return q;
        };

        std::function<bool(Point<DRootTwo>)> test = [=](Point<DRootTwo> p)
        {
            DRootTwo x, y;
            std::tie(x, y) = p;
            return x * x + y * y <= s;
        };

        return ConvexSet<T>(el, test, intersect);
    }

    template <typename T>
    Operator<T> op_fromDRootTwo(Operator<DRootTwo> op)
    {
        Operator<T> result;
        for (size_t i = 0; i < op.size1(); i++)
        {
            for (size_t j = 0; j < op.size2(); j++)
            {
                result(i, j) = ring::fromDRootTwo<T>(op(i, j));
            }
        }
        return result;
    }

    template <typename T>
    std::tuple<T, double> operator_to_bz(Operator<T> op)
    {
        T a = op(0, 0);
        T b = op(0, 1);
        T d = op(1, 1);
        T lz = d * ring::recip<T>(a);
        double z = 0.5 * logBase_double<T>(lambda<T>(), lz);
        return std::make_tuple(b, z);
    }

    template <typename T>
    Pair<T> operator_to_bl2z(Operator<T> op)
    {
        T a = op(0, 0);
        T b = op(0, 1);
        T d = op(1, 1);
        T l2z = d * ring::recip<T>(a);
        return std::make_tuple(b, l2z);
    }

    template <typename T>
    T det(Operator<T> op)
    {
        return op(0, 0) * op(1, 1) - op(0, 1) * op(1, 0);
    }

    template <typename T>
    T operator_skew(Operator<T> op)
    {
        return op(0, 1) * op(1, 0);
    }

    template <typename T>
    T uprightness(Operator<T> op)
    {
        return constants::pi<T>() / 4 * bmp::sqrt(det(op) / (op(0, 0) * op(1, 1)));
    }

    template <typename T>
    T skew(OperatorPair<T> pair)
    {
        Operator<T> op1 = std::get<0>(pair);
        Operator<T> op2 = std::get<1>(pair);
        return operator_skew(op1) + operator_skew(op2);
    }

    template <typename T>
    double bias(OperatorPair<T> pair)
    {
        Operator<T> op1 = std::get<0>(pair);
        Operator<T> op2 = std::get<1>(pair);
        double z = std::get<1>(operator_to_bz(op1));
        double zeta = std::get<1>(operator_to_bz(op2));
        return zeta - z;
    }

    template <typename T>
    Operator<T> opR() { return ring::roothalf<T>() * mat::matrix2x2<T>(1, -1, 1, 1); }

    template <typename T>
    Operator<T> opA() { return mat::matrix2x2<T>(1, -2, 0, 1); }

    template <typename T>
    Operator<T> opA_inv() { return mat::matrix2x2<T>(1, 2, 0, 1); }

    template <typename T>
    Operator<T> opA_power(Integer k)
    {
        return (k >= 0) ? operator_pow<T>(opA<T>(), k) : operator_pow<T>(opA_inv<T>(), -k);
    }

    template <typename T>
    Operator<T> opB() { return mat::matrix2x2<T>(1, ring::roottwo<T>(), 0, 1); }

    template <typename T>
    Operator<T> opB_inv() { return mat::matrix2x2<T>(1, -ring::roottwo<T>(), 0, 1); }

    template <typename T>
    Operator<T> opB_power(Integer k)
    {
        return (k >= 0) ? operator_pow<T>(opB<T>(), k) : operator_pow<T>(opB_inv<T>(), -k);
    }

    template <typename T>
    Operator<T> opK() { return ring::roothalf<T>() * mat::matrix2x2<T>(-lambda_inv<T>(), -1, lambda<T>(), 1); }

    template <typename T>
    Operator<T> opX() { return mat::matrix2x2<T>(0, 1, 1, 0); }

    template <typename T>
    Operator<T> opZ() { return mat::matrix2x2<T>(1, 0, 0, -1); }

    template <typename T>
    Operator<T> opS() { return mat::matrix2x2<T>(lambda<T>(), 0, 0, lambda_inv<T>()); }

    template <typename T>
    Operator<T> opS_inv() { return mat::matrix2x2<T>(lambda_inv<T>(), 0, 0, lambda<T>()); }

    template <typename T>
    Operator<T> opS_power(Integer k)
    {
        return (k >= 0) ? operator_pow<T>(opS<T>(), k) : operator_pow<T>(opS_inv<T>(), -k);
    }

    template <typename T>
    OperatorPair<T> action(OperatorPair<T> pair, Operator<DRootTwo> g)
    {
        Operator<T> a = std::get<0>(pair);
        Operator<T> b = std::get<1>(pair);
        Operator<T> g2 = op_fromDRootTwo<T>(g);
        Operator<T> g4 = op_fromDRootTwo<T>(adj2(g));
        Operator<T> g1 = adj(g2);
        Operator<T> g3 = adj(g4);
        Operator<T> p1 = prod(g1, a);
        Operator<T> m1 = prod(p1, g2);
        Operator<T> p2 = prod(g3, b);
        Operator<T> m2 = prod(p2, g4);
        return std::make_tuple(m1, m2);
    }

    template <typename T>
    Operator<T> shift_sigma(Integer k, Operator<T> m)
    {
        Operator<T> op;
        op(0, 0) = lambdapower<T>(k) * m(0, 0);
        op(0, 1) = m(0, 1);
        op(1, 0) = m(1, 0);
        op(1, 1) = lambdapower<T>(-k) * m(1, 1);
        return op;
    }

    template <typename T>
    Operator<T> shift_tau(Integer k, Operator<T> m)
    {
        Operator<T> op;
        op(0, 0) = lambdapower<T>(-k) * m(0, 0);
        op(0, 1) = signpower<T>(k) * m(0, 1);
        op(1, 0) = m(1, 0) * signpower<T>(k);
        op(1, 1) = lambdapower<T>(k) * m(1, 1);
        return op;
    }

    template <typename T>
    OperatorPair<T> shift_state(Integer k, OperatorPair<T> pair)
    {
        Operator<T> d = std::get<0>(pair);
        Operator<T> delta = std::get<1>(pair);
        return std::make_tuple(shift_sigma(k, d), shift_tau(k, delta));
    }

    template <typename T>
    Integer lemma_A(T z, T zeta)
    {
        T c = std::min<T>(z, zeta);
        return std::max<Integer>(1, ring::floor_of(bmp::pow(lambda<T>(), c) / T(2)));
    }

    template <typename T>
    Integer lemma_B(T z, T zeta)
    {
        T c = std::min<T>(z, zeta);
        return std::max<Integer>(1, ring::floor_of(bmp::pow(lambda<T>(), c) / ring::roottwo<T>()));
    }

    template <typename T>
    Integer lemma_A_l2(T l2z, T l2zeta)
    {
        T l2c = std::min<T>(l2z, l2zeta);
        return std::max<Integer>(1, ring::intsqrt(ring::floor_of(l2c / 4)));
    }

    template <typename T>
    Integer lemma_B_l2(T l2z, T l2zeta)
    {
        T l2c = std::min<T>(l2z, l2zeta);
        return std::max<Integer>(1, ring::intsqrt(ring::floor_of(l2c / 2)));
    }

    template <typename T>
    Maybe<Operator<DRootTwo>> step_lemma(OperatorPair<T> pair)
    {
        Operator<T> matA = std::get<0>(pair);
        Operator<T> matB = std::get<1>(pair);
        T b, l2z;
        std::tie(b, l2z) = operator_to_bl2z(matA);
        T beta, l2zeta;
        std::tie(beta, l2zeta) = operator_to_bl2z(matB);

        auto logLambda = [](T a)
        { return logBase_double<T>(lambda<T>(), a); };

        T l2z_minus_zeta = l2z / l2zeta;

        auto wlog_using = [=](Operator<DRootTwo> op) -> Maybe<Operator<DRootTwo>>
        {
            Operator<T> matA2, matB2;
            std::tie(matA2, matB2) = action(std::make_tuple(matA, matB), op);
            Maybe<Operator<DRootTwo>> maybe_op2 = step_lemma(std::make_tuple(matA2, matB2));
            Maybe<Operator<DRootTwo>> result;
            if (!maybe_op2.has_value())
            {
                result = op;
            }
            else
            {
                result = prod(op, maybe_op2.value());
            }
            return result;
        };

        auto with_shift = [=](Integer k) -> Maybe<Operator<DRootTwo>>
        {
            Operator<T> matA2, matB2;
            std::tie(matA2, matB2) = shift_state(k, std::make_tuple(matA, matB));
            Maybe<Operator<DRootTwo>> maybe_op2 = step_lemma(std::make_tuple(matA2, matB2));
            Maybe<Operator<DRootTwo>> result;
            if (!maybe_op2.has_value())
            {
                result = std::nullopt;
            }
            else
            {
                result = shift_sigma(k, maybe_op2.value());
            }
            return result;
        };

        if (beta < 0)
        {
            return wlog_using(opZ<DRootTwo>());
        }
        if (l2z * l2zeta < 1)
        {
            return wlog_using(opX<DRootTwo>());
        }
        if (l2z_minus_zeta > T("33.971") || l2z_minus_zeta < T("0.029437"))
        {
            return wlog_using(opS_power<DRootTwo>(std::round(logLambda(l2z_minus_zeta) / 8)));
        }
        if (skew(std::make_tuple(matA, matB)) <= 15)
        {
            return std::nullopt;
        }
        if (l2z_minus_zeta > T("5.8285") || l2z_minus_zeta < T("0.17157"))
        {
            return with_shift(round(logLambda(l2z_minus_zeta) / 4));
        }
        if (within<T>(l2z, T("0.24410"), T("4.0968")) && within<T>(l2zeta, T("0.24410"), T("4.0968")))
        {
            return opR<DRootTwo>();
        }
        if (b >= 0 && l2z <= T("1.6969"))
        {
            return opK<DRootTwo>();
        }
        if (b >= 0 && l2zeta <= T("1.6969"))
        {
            return adj2<DRootTwo>(opK<DRootTwo>());
        }
        if (b >= 0)
        {
            return opA_power<DRootTwo>(lemma_A_l2<T>(l2z, l2zeta));
        }
        return opB_power<DRootTwo>(lemma_B_l2<T>(l2z, l2zeta));
    }

    template <typename T>
    Operator<DRootTwo> reduction(OperatorPair<T> st)
    {
        Maybe<Operator<DRootTwo>> sl = step_lemma(st);
        if (!sl.has_value())
        {
            return mat::matrix2x2<DRootTwo>(DRootTwo(1), DRootTwo(0), DRootTwo(0), DRootTwo(1));
        }
        Operator<DRootTwo> opG = sl.value();
        Operator<DRootTwo> opG2 = reduction<T>(action<T>(st, opG));
        return prod(opG, opG2);
    }

    template <typename T>
    Operator<DRootTwo> to_upright(OperatorPair<T> pair)
    {
        Operator<T> a = std::get<0>(pair);
        Operator<T> b = std::get<1>(pair);
        Operator<T> a2 = a / bmp::sqrt(det<T>(a));
        Operator<T> b2 = b / bmp::sqrt(det<T>(b));
        return reduction<T>(std::make_tuple(a2, b2));
    }

    template <typename T>
    Operator<DRootTwo> to_upright_sets(ConvexSet<T> setA, ConvexSet<T> setB)
    {
        return to_upright(std::make_tuple(setA.el().op(), setB.el().op()));
    }

    template <typename T>
    Point<T> point_transform(Operator<T> opG, Point<T> p)
    {
        T x, y;
        x = std::get<0>(p);
        y = std::get<1>(p);
        T a, b, c, d;
        a = opG(0, 0);
        b = opG(0, 1);
        c = opG(1, 0);
        d = opG(1, 1);
        return std::make_tuple(a * x + b * y, c * x + d * y);
    }

    template <typename T>
    Ellipse<T> ellipse_transform(Operator<T> opG, Ellipse<T> ell)
    {
        Operator<T> opG_inv = special_inverse<T>(opG);
        Operator<T> p1 = prod(adj(opG_inv), ell.op());
        Operator<T> new_op = prod(p1, opG_inv);
        Point<T> new_center = point_transform(opG, ell.p());
        return Ellipse<T>(new_op, new_center);
    }

    CharFun charfun_transform(Operator<DRootTwo> opG, CharFun f)
    {
        Operator<DRootTwo> opG_inv = special_inverse<DRootTwo>(opG);
        return [=](Point<DRootTwo> p)
        {
            return f(point_transform(opG_inv, p));
        };
    }

    template <typename T>
    LineIntersector<T> lineintersector_transform(Operator<DRootTwo> opG, LineIntersector<T> intA)
    {
        Operator<DRootTwo> opG_inv = special_inverse<DRootTwo>(opG);
        return [=](Point<DRootTwo> v2, Point<DRootTwo> w2)
        {
            Point<DRootTwo> v = point_transform(opG_inv, v2);
            Point<DRootTwo> w = point_transform(opG_inv, w2);
            return intA(v, w);
        };
    }

    template <typename T>
    ConvexSet<T> convex_transform(Operator<DRootTwo> opG, ConvexSet<T> set)
    {
        Ellipse<T> new_ell = ellipse_transform<T>(op_fromDRootTwo<T>(opG), set.el());
        LineIntersector<T> new_int = lineintersector_transform<T>(opG, set.intersect_fun());
        CharFun new_test = charfun_transform(opG, set.test_fun());
        return ConvexSet<T>(new_ell, new_test, new_int);
    }

    template <typename T>
    Tuple2By2<T> boundingbox_ellipse(Ellipse<T> el)
    {
        T x, y;
        Point<T> ctrA = el.p();
        x = std::get<0>(ctrA);
        y = std::get<1>(ctrA);
        Operator<T> matA = el.op();
        T a, b, c, d;
        a = matA(0, 0);
        b = matA(0, 1);
        c = matA(1, 0);
        d = matA(1, 1);
        T sqrt_det = bmp::sqrt(det(matA));
        T w = bmp::sqrt(d) / sqrt_det;
        T h = bmp::sqrt(a) / sqrt_det;
        return std::make_tuple(std::make_tuple(x - w, x + w), std::make_tuple(y - h, y + h));
    }

    template <typename T>
    Tuple2By2<T> boundingbox(ConvexSet<T> set)
    {
        return boundingbox_ellipse(set.el());
    }

    template <typename T>
    std::function<List<DOmega>(Integer)> gridpoints2_scaled_with_gridop(
        ConvexSet<T> setA, ConvexSet<T> setB, Operator<DRootTwo> opG)
    {
        Operator<DRootTwo> opG_inv = special_inverse<DRootTwo>(opG);
        ConvexSet<T> setA_prime = convex_transform(opG_inv, setA);
        ConvexSet<T> setB_prime = convex_transform(adj2<DRootTwo>(opG_inv), setB);
        Tuple2By2<T> bboxA_prime = boundingbox(setA_prime);
        Tuple2By2<T> bboxB_prime = boundingbox(setB_prime);
        T x0A, x1A, y0A, y1A;
        std::tie(x0A, x1A) = std::get<0>(bboxA_prime);
        std::tie(y0A, y1A) = std::get<1>(bboxA_prime);
        T x0B, x1B, y0B, y1B;
        std::tie(x0B, x1B) = std::get<0>(bboxB_prime);
        std::tie(y0B, y1B) = std::get<1>(bboxB_prime);

        std::function<List<DOmega>(Integer)> solutions_fun = [=](Integer k)
        {
            if (k < 0)
            {
                throw std::invalid_argument("k >= 0 is required");
            }

            Pair<T> intervalA = std::make_tuple(y0A, y1A);
            Pair<T> intervalB = std::make_tuple(y0B, y1B);

            List<DRootTwo> xs = gridpoints_scaled<T>(
                x0A, x1A + lambda<T>(), x0B, x1B + lambda<T>(), k + 1);

            List<DOmega> result(0);

            if (xs.size() == 0)
            {
                return result;
            }

            DRootTwo x0 = xs.at(0);
            DRootTwo x0_bul = x0.adj2();
            DRootTwo dx = ring::pow_non_neg(ring::roothalf<DRootTwo>(), k);
            DRootTwo dx_bul = ring::adj2(dx);

            List<DRootTwo> beta_prime_list = gridpoints_scaled(
                fatten_interval(intervalA), fatten_interval(intervalB), k + 1);

            for (DRootTwo beta_prime : beta_prime_list)
            {
                DRootTwo beta_prime_bul = ring::adj2<DRootTwo>(beta_prime);

                Point<DRootTwo> p1A = std::make_tuple(x0, beta_prime);
                Point<DRootTwo> p2A = std::make_tuple(dx, DRootTwo(0));
                Maybe<Pair<T>> iA = setA_prime.intersect(p1A, p2A);

                Point<DRootTwo> p1B = std::make_tuple(x0_bul, beta_prime_bul);
                Point<DRootTwo> p2B = std::make_tuple(dx_bul, DRootTwo(0));
                Maybe<Pair<T>> iB = setB_prime.intersect(p1B, p2B);

                if (!iA.has_value() || !iB.has_value())
                {
                    continue;
                }

                T t0A, t1A, t0B, t1B;
                std::tie(t0A, t1A) = iA.value();
                std::tie(t0B, t1B) = iB.value();

                T dtA = T(10) * ring::recip<T>(std::max<T>(
                                    T(10), ring::pow_non_neg<T>(T(2), k) * (t1B - t0B)));
                T dtB = T(10) * ring::recip<T>(std::max<T>(
                                    T(10), ring::pow_non_neg<T>(T(2), k) * (t1A - t0A)));

                DRootTwo rk = ring::pow_non_neg<DRootTwo>(ring::roottwo<DRootTwo>(), k);
                List<DRootTwo> alpha_prime_offs_list = gridpoints_scaled_parity(
                    (beta_prime - x0) * rk, t0A - dtA, t1A + dtA, t0B - dtB, t1B + dtB, 1);

                for (DRootTwo alpha_prime_offs : alpha_prime_offs_list)
                {
                    DRootTwo alpha_prime = alpha_prime_offs * dx + x0;
                    DRootTwo alpha, beta;
                    std::tie(alpha, beta) = point_transform(opG, std::make_tuple(alpha_prime, beta_prime));
                    Point<DRootTwo> p = std::make_tuple(alpha, beta);
                    Point<DRootTwo> p2 = std::make_tuple(alpha.adj2(), beta.adj2());
                    if (setA.test(p) && setB.test(p2))
                    {
                        DOmega i = ring::i<DOmega>();
                        DOmega z = ring::fromDRootTwo<DOmega>(alpha) + i * ring::fromDRootTwo<DOmega>(beta);
                        result.push_back(z);
                    }
                }
            }
            return result;
        };
        return solutions_fun;
    }

    template <typename T>
    std::function<List<DOmega>(Integer)> gridpoints2_scaled(ConvexSet<T> setA, ConvexSet<T> setB)
    {
        Operator<DRootTwo> opG = to_upright_sets<T>(setA, setB);
        return gridpoints2_scaled_with_gridop(setA, setB, opG);
    }

    /**
     * @brief Modified version of the Haskell implementation.
     *
     * The Haskell implementation returns an infinite list of tuples of the form (k, l_k),
     * where k is a non-negative integer and l_k is a list of the solutions to the scaled
     * grid problem that first appear for a scale of k. This kind of infinite list
     * wouldn't work directly in C++ because C++ doesn't use lazy evaluation. Instead, this
     * function returns a new function that takes in k and returns the list of solutions
     * that first appear for that k.
     */
    template <typename T>
    std::function<List<DOmega>(Integer)> gridpoints2_increasing_with_gridop(
        ConvexSet<T> setA, ConvexSet<T> setB, Operator<DRootTwo> opG)
    {
        std::function<List<DOmega>(Integer)> solutions_fun = gridpoints2_scaled_with_gridop(setA, setB, opG);
        std::function<List<DOmega>(Integer)> exact_solutions = [=](Integer k)
        {
            if (k == 0)
            {
                return solutions_fun(0);
            }
            List<DOmega> all_sols = solutions_fun(k);
            List<DOmega> exact;
            std::copy_if(all_sols.begin(), all_sols.end(), std::back_inserter(exact), [=](DOmega d)
                         { return ring::denomexp(d) == k; });
            return exact;
        };
        return exact_solutions;
    }

    /**
     * @brief Modified version of the Haskell implementation.
     *
     * See the documentation for gridpoints2_increasing_with_gridop: the same explanation
     * about the difference in return type for Haskell vs. C++ applies here.
     */
    template <typename T>
    std::function<List<DOmega>(Integer)> gridpoints2_increasing(ConvexSet<T> setA, ConvexSet<T> setB)
    {
        return gridpoints2_increasing_with_gridop(setA, setB, to_upright_sets(setA, setB));
    }
}