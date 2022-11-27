#include "ring.h"
#include "quadratic.h"
#include <vector>

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
CharFun ConvexSet<T>::test() const
{
    return test_;
}

template <typename T>
LineIntersector<T> ConvexSet<T>::intersect() const
{
    return intersect_;
}

namespace gridprob
{
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
}