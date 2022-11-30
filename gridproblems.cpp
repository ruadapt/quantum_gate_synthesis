#include "ring.h"
#include "quadratic.h"
#include <vector>
#include <cmath>

template <typename T>
T fst(Point<T> p)
{
    return std::get<0>(p);
}

template <typename T>
T snd(Point<T> p)
{
    return std::get<1>(p);
}

template <typename T>
Operator<T> Ellipse<T>::op() const
{
    return op_;
}

template <typename T>
Point<T> Ellipse<T>::p() const
{
    return p_;
}

template <typename T>
Ellipse<T> ConvexSet<T>::el() const
{
    return el_;
}

template <typename T>
bool ConvexSet<T>::test(Point<DRootTwo> p) const
{
    return test_(p);
}

template <typename T>
std::optional<std::tuple<T, T>> ConvexSet<T>::intersect(Point<DRootTwo> p1, Point<DRootTwo> p2) const
{
    return intersect_(p1, p2);
}

namespace gridprob
{
    template <typename T, int M, int N>
    Matrix<T, M, N> matrixPow(Matrix<T, M, N> base, int exp)
    {
        assert(exp >= 0);
        Matrix<T, M, N> result;
        if (exp == 0)
        {
            result = 1 + result;
            return result;
        }
        if (exp == 1)
        {
            return base;
        }
        if (exp % 2 == 0)
        {
            T newBase = base * base;
            int newExp = exp / 2;
            return matrixPow(newBase, newExp);
        }
        T newBase = base * base;
        int newExp = (exp - 1) / 2;
        return prod(base, matrixPow(newBase, newExp));
    }

    template <typename T>
    Operator<T> operatorPow(Operator<T> base, int exp)
    {
        return matrixPow<T, 2, 2>(base, exp);
    }

    double logBase(double b, double x)
    {
        return log(x) / log(b);
    }

    template <typename T>
    T lambda()
    {
        return T(1) + ring::rootTwo<T>();
    }

    template <typename T>
    T lambdaInv()
    {
        return ring::rootTwo<T>() - 1;
    }

    template <typename T>
    bool within(const T x, const T low, const T high)
    {
        return (low <= x && x <= high);
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
    double logBaseDouble(T b, T x)
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
        return logBaseDouble(ring::recip(b), x);
    }

    template <typename T>
    T iprod(Point<T> p1, Point<T> p2)
    {
        return fst(p1) * fst(p2) + snd(p1) * snd(p2);
    }

    template <typename T>
    std::vector<ZRootTwo> gridpointsInternal(T x0, T x1, T y0, T y1)
    {
        T dx = x1 - x0;
        T dy = y1 - y0;

        if (dy < 0 || dx < 0)
        {
            return std::vector<ZRootTwo>{};
        }

        auto baseCase = [=]
        {
            Integer amin = ring::ceiling_of((x0 + y0) * ring::recip<T>(2));
            Integer amax = ring::floor_of((x1 + y1) * ring::recip<T>(2));
            std::vector<ZRootTwo> results;
            for (Integer a = amin; a <= amax; a++)
            {
                T aT = ring::fromInteger<T>(a);
                Integer bmin = ring::ceiling_of((aT - y1) * ring::recip(ring::rootTwo<T>()));
                Integer bmax = ring::floor_of((aT - y0) * ring::recip(ring::rootTwo<T>()));
                for (Integer b = bmin; b <= bmax; b++)
                {
                    results.push_back(ZRootTwo(a, b));
                }
            }
            return results;
        };

        if (dy == 0 && dx == 0)
        {
            return baseCase();
        }

        if (dy == 0 && dx > 0)
        {
            std::vector<ZRootTwo> subResults = gridpointsInternal(y0, y1, x0, x1);
            std::vector<ZRootTwo> results(subResults.size());
            std::transform(subResults.begin(), subResults.end(), results.begin(), [](ZRootTwo z)
                           { return z.adj2(); });
            return results;
        }

        // floorlog requires a positive argument, so we don't compute it until
        // we've ruled out the possibility that dy <= 0.
        int n = ring::mpzToInt(std::get<0>(floorlog<T>(lambda<T>(), dy)));
        int m = -n;

        T lT = lambda<T>();
        T lInvT = lambdaInv<T>();
        T lambdaMT = ring::powInt(lT, m);
        T lambdaNT = ring::powInt(lT, n);
        T lambdaBulNT = ring::powInt(-lInvT, n);
        T lambdaInvMT = ring::powInt(lInvT, m);
        T lambdaBulInvMT = ring::powInt(-lT, m);
        T lambdaInvNT = ring::powInt(lInvT, n);
        ZRootTwo lZ = lambda<ZRootTwo>();
        ZRootTwo lInvZ = lambdaInv<ZRootTwo>();

        if (dy >= lT && ring::even(n))
        {
            std::vector<ZRootTwo> subResults = gridpointsInternal(
                lambdaNT * x0, lambdaNT * x1, lambdaBulNT * y0, lambdaBulNT * y1);
            std::vector<ZRootTwo> results(subResults.size());
            ZRootTwo lambdaInvNZ = ring::powNonNeg(lInvZ, n);
            std::transform(subResults.begin(), subResults.end(), results.begin(),
                           [&lambdaInvNZ = std::as_const(lambdaInvNZ)](ZRootTwo z)
                           { return lambdaInvNZ * z; });
            return results;
        }

        if (dy >= lT && ring::odd(n))
        {
            std::vector<ZRootTwo> subResults = gridpointsInternal(
                lambdaNT * x0, lambdaNT * x1, lambdaBulNT * y1, lambdaBulNT * y0);
            std::vector<ZRootTwo> results(subResults.size());
            ZRootTwo lambdaInvNZ = ring::powNonNeg(lInvZ, n);
            std::transform(subResults.begin(), subResults.end(), results.begin(),
                           [&lambdaInvNZ = std::as_const(lambdaInvNZ)](ZRootTwo z)
                           { return lambdaInvNZ * z; });
            return results;
        }

        if (dy > 0 && dy < 1 && ring::even(n))
        {
            std::vector<ZRootTwo> subResults = gridpointsInternal(
                lambdaInvMT * x0, lambdaInvMT * x1, lambdaBulInvMT * y0, lambdaBulInvMT * y1);
            std::vector<ZRootTwo> results(subResults.size());
            ZRootTwo lambdaMZ = ring::powNonNeg(lZ, m);
            std::transform(subResults.begin(), subResults.end(), results.begin(),
                           [&lambdaMZ = std::as_const(lambdaMZ)](ZRootTwo z)
                           { return lambdaMZ * z; });
            return results;
        }

        if (dy > 0 && dy < 1 && ring::odd(n))
        {
            std::vector<ZRootTwo> subResults = gridpointsInternal(
                lambdaInvMT * x0, lambdaInvMT * x1, lambdaBulInvMT * y1, lambdaBulInvMT * y0);
            std::vector<ZRootTwo> results(subResults.size());
            ZRootTwo lambdaMZ = ring::powNonNeg(lZ, m);
            std::transform(subResults.begin(), subResults.end(), results.begin(),
                           [&lambdaMZ = std::as_const(lambdaMZ)](ZRootTwo z)
                           { return lambdaMZ * z; });
            return results;
        }

        return baseCase();
    }

    template <typename T>
    std::vector<ZRootTwo> gridpoints(T x0, T x1, T y0, T y1)
    {
        Integer a = ring::div(ring::floor_of(x0 + y0), 2_mpz);
        Integer b = ring::div(ring::floor_of(ring::rootTwo<T>() * (x0 - y0)), 4_mpz);
        ZRootTwo alpha = ZRootTwo(a, b);
        T xoff = ring::fromZRootTwo<T>(alpha);
        T yoff = ring::fromZRootTwo<T>(ring::adj2(alpha));
        T x0new = x0 - xoff;
        T x1new = x1 - xoff;
        T y0new = y0 - yoff;
        T y1new = y1 - yoff;

        std::vector<ZRootTwo> rawPoints = gridpointsInternal(x0new, x1new, y0new, y1new);
        std::vector<ZRootTwo> shifted(rawPoints.size());
        std::transform(rawPoints.begin(), rawPoints.end(), shifted.begin(), [=](ZRootTwo z)
                       { return z + alpha; });
        std::vector<ZRootTwo> results;
        std::copy_if(shifted.begin(), shifted.end(), std::back_inserter(results), [=](ZRootTwo z)
                     { 
                        T zz = ring::fromZRootTwo<T>(z);
                        return within<T>(zz, x0, x1) && within<T>(ring::adj2(zz), y0, y1); });
        return results;
    }

    template <typename T>
    std::vector<DRootTwo> gridpointsScaled(T x0, T x1, T y0, T y1, Integer k)
    {
        if (k < 0)
        {
            throw std::invalid_argument("gridpointsScaled: k >= 0 is required");
        }
        T scaleT = ring::powNonNeg(ring::rootHalf<T>(), k);
        T scaleInvT = ring::powNonNeg(ring::rootTwo<T>(), k);
        DRootTwo scaleD = ring::powNonNeg(ring::rootHalf<DRootTwo>(), k);
        DRootTwo scaleInvD = ring::powNonNeg(ring::rootTwo<DRootTwo>(), k);
        T x0new = scaleInvT * x0;
        T x1new = scaleInvT * x1;
        T y0new = ring::even(k) ? (scaleInvT * y0) : (-scaleInvT * y1);
        T y1new = ring::even(k) ? (scaleInvT * y1) : (-scaleInvT * y0);
        std::vector<ZRootTwo> w = gridpoints(x0new, x1new, y0new, y1new);
        std::vector<DRootTwo> results(w.size());
        std::transform(w.begin(), w.end(), results.begin(), [=](ZRootTwo z)
                       { return scaleD * ring::fromZRootTwo<DRootTwo>(z); });
        return results;
    }

    template <typename T>
    std::vector<DRootTwo> gridpointsScaledParity(DRootTwo beta, T x0, T x1, T y0, T y1, Integer k)
    {
        if (k < 1)
        {
            throw std::invalid_argument("gridpointsScaledParity: k >= 1 is required");
        }
        if (ring::denomExp(beta) <= k - 1)
        {
            return gridpointsScaled(x0, x1, y0, y1, k - 1);
        }
        DRootTwo offs = ring::powNonNeg(ring::rootHalf<DRootTwo>(), k);
        DRootTwo offsBul = ring::adj2(offs);
        T offsPrime = ring::fromDRootTwo<T>(offs);
        T offsBulPrime = ring::fromDRootTwo<T>(offsBul);
        std::vector<DRootTwo> z = gridpointsScaled(x0 + offsPrime, x1 + offsPrime, y0 + offsBulPrime, y1 + offsBulPrime, k - 1);
        std::vector<DRootTwo> results(z.size());
        std::transform(z.begin(), z.end(), results.begin(), [=](DRootTwo d)
                       { return d - offs; });
        return results;
    }

    template <typename T>
    Point<T> pointFromDRootTwo(Point<DRootTwo> p)
    {
        return std::make_tuple(ring::fromDRootTwo<T>(fst(p)), ring::fromDRootTwo<T>(snd(p)));
    }

    template <typename T>
    Operator<T> makeOperator(T x0, T x1, T x2, T x3)
    {
        Operator<T> op;
        op(0, 0) = x0;
        op(0, 1) = x1;
        op(1, 0) = x2;
        op(1, 1) = x3;
        return op;
    }

    template <typename T>
    ConvexSet<T> unitDisk()
    {
        Operator<T> op = makeOperator<T>(1, 0, 0, 1);
        Point<T> p = std::make_tuple<T, T>(0, 0);
        Ellipse<T> el = Ellipse<T>(op, p);

        std::function<std::optional<std::tuple<T, T>>(Point<DRootTwo>, Point<DRootTwo>)> intersect = [](Point<DRootTwo> p, Point<DRootTwo> v)
        {
            QRootTwo a = ring::fromDRootTwo<QRootTwo>(iprod<DRootTwo>(v, v));
            QRootTwo b = ring::fromDRootTwo<QRootTwo>(DRootTwo(2) * iprod<DRootTwo>(v, p));
            QRootTwo c = ring::fromDRootTwo<QRootTwo>(iprod<DRootTwo>(p, p) - 1);
            std::optional<std::tuple<Real, Real>> q = quadratic<QRootTwo>(a, b, c);
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
        Operator<T> op = makeOperator<T>(r, 0, 0, r);
        Point<T> p = std::make_tuple<T, T>(0, 0);
        Ellipse<T> el = Ellipse<T>(op, p);

        std::function<std::optional<std::tuple<T, T>>(Point<DRootTwo>, Point<DRootTwo>)> intersect = [=](Point<DRootTwo> p, Point<DRootTwo> v)
        {
            QRootTwo a = ring::fromDRootTwo<QRootTwo>(iprod<DRootTwo>(v, v));
            QRootTwo b = ring::fromDRootTwo<QRootTwo>(DRootTwo(2) * iprod<DRootTwo>(v, p));
            QRootTwo c = ring::fromDRootTwo<QRootTwo>(iprod<DRootTwo>(p, p) - s);
            std::optional<std::tuple<Real, Real>> q = quadratic<QRootTwo>(a, b, c);
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
    Operator<T> opFromDRootTwo(Operator<DRootTwo> op)
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
    std::tuple<T, double> operatorToBz(Operator<T> op)
    {
        T a = op(0, 0);
        T b = op(0, 1);
        T d = op(1, 1);
        T lz = d * ring::recip<T>(a);
        return 0.5 * logBaseDouble(lambda<T>(), lz);
    }

    template <typename T>
    T det(Operator<T> op)
    {
        return op(0, 0) * op(1, 1) - op(0, 1) * op(1, 0);
    }

    template <typename T>
    T operatorSkew(Operator<T> op)
    {
        return op(0, 1) * op(1, 0);
    }

    template <typename T>
    T uprightness(Operator<T> op)
    {
        return M_PI / 4 * sqrt(det(op) / (op(0, 0) * op(1, 1)));
    }

    template <typename T>
    T skew(OperatorPair<T> pair)
    {
        Operator<T> op1 = std::get<0>(pair);
        Operator<T> op2 = std::get<1>(pair);
        return operatorSkew(op1) + operatorSkew(op2);
    }

    template <typename T>
    double bias(OperatorPair<T> pair)
    {
        Operator<T> op1 = std::get<0>(pair);
        Operator<T> op2 = std::get<1>(pair);
        double z = std::get<1>(operatorToBz(op1));
        double zeta = std::get<1>(operatorToBz(op2));
        return zeta - z;
    }

    template <typename T>
    Operator<T> opR() { return ring::rootHalf<T>() * makeOperator<T>(1, -1, 1, 1); }

    template <typename T>
    Operator<T> opA() { return makeOperator<T>(1, -2, 0, 1); }

    template <typename T>
    Operator<T> opAInv() { return makeOperator<T>(1, 2, 0, 1); }

    template <typename T>
    Operator<T> opAPower(Integer k)
    {
        return (k >= 0) ? operatorPow<T>(opA<T>(), k) : operatorPow<T>(opAInv<T>(), -k);
    }

    template <typename T>
    Operator<T> opB() { return makeOperator<T>(1, ring::rootTwo<T>(), 0, 1); }

    template <typename T>
    Operator<T> opBInv() { return makeOperator<T>(1, -ring::rootTwo<T>(), 0, 1); }

    template <typename T>
    Operator<T> opBPower(Integer k)
    {
        return (k >= 0) ? operatorPow<T>(opB<T>(), k) : operatorPow<T>(opBInv<T>(), -k);
    }

    template <typename T>
    Operator<T> opK() { return ring::rootHalf<T>() * makeOperator<T>(-lambdaInv<T>(), -1, lambda<T>(), 1); }

    template <typename T>
    Operator<T> opX() { return makeOperator<T>(0, 1, 1, 0); }

    template <typename T>
    Operator<T> opZ() { return makeOperator<T>(1, 0, 0, -1); }

    template <typename T>
    Operator<T> opS() { return makeOperator<T>(lambda<T>(), 0, 0, lambdaInv<T>()); }

    template <typename T>
    Operator<T> opSInv() { return makeOperator<T>(lambdaInv<T>(), 0, 0, lambda<T>()); }

    template <typename T>
    Operator<T> opSPower(Integer k)
    {
        return (k >= 0) ? operatorPow<T>(opS<T>(), k) : operatorPow<T>(opSInv<T>(), -k);
    }
};