#include "gridSynth.h"
#include "quadratic.h"
#include "ring.h"

namespace gridsynth
{
    namespace gp = gridprob;
    namespace bmp = boost::multiprecision;

    template <typename T>
    ConvexSet<T> epsilon_region(T epsilon, T theta)
    {
        T d = 1 - bmp::pow(epsilon, 2) / 2;
        T zx = bmp::cos(-theta / 2);
        T zy = bmp::sin(-theta / 2);
        Point<T> z = Point<T>{zx, zy};
        T ev1 = 4 * bmp::pow(1 / epsilon, 4);
        T ev2 = bmp::pow(1 / epsilon, 2);
        Operator<T> mmat = matrix2x2<T>(ev1, 0, 0, ev2);
        Operator<T> bmat = matrix2x2<T>(zx, -zy, zy, zx);
        Operator<T> temp1 = prod(bmat, mmat);
        Operator<T> mat = prod(temp1, gp::special_inverse(bmat));
        Point<T> ctr = Point<T>{d * zx, d * zy};
        Ellipse<T> ell = Ellipse<T>(mat, ctr);

        auto intersect = [=](Point<DRootTwo> p, Point<DRootTwo> v) -> Maybe<Pair<T>>
        {
            DRootTwo a = gp::iprod(v, v);
            DRootTwo b = DRootTwo(2) * gp::iprod(v, p);
            DRootTwo c = gp::iprod(p, p) - 1;
            Maybe<Pair<T>> q = quadratic(ring::fromDRootTwo<QRootTwo>(a), ring::fromDRootTwo<QRootTwo>(b), ring::fromDRootTwo<QRootTwo>(c));

            if (!q.has_value())
            {
                return Maybe<Pair<T>>();
            }
            T t0, t1;
            std::tie(t0, t1) = q.value();
            T vz = gp::iprod(gp::pointFromDRootTwo<T>(v), z);
            T rhs = d - gp::iprod(gp::pointFromDRootTwo<T>(p), z);
            if (vz == 0 && rhs <= 0)
            {
                return Maybe<Pair<T>>(Pair<T>{t0, t1});
            }
            if (vz == 0)
            {
                return Maybe<Pair<T>>();
            }
            T t2 = rhs / vz;
            if (vz > 0)
            {
                return Maybe<Pair<T>>(Pair<T>{std::max(t0, t2), t1});
            }
            return Maybe<Pair<T>>(Pair<T>{t0, std::min(t1, t2)});
        };

        auto tst = [=](Point<DRootTwo> p) -> bool
        {
            DRootTwo x, y;
            std::tie(x, y) = p;
            bool cond1 = x * x + y * y <= DRootTwo(1);
            bool cond2 = zx * ring::fromDRootTwo<T>(x) + zy * ring::fromDRootTwo<T>(y) >= d;
            return cond1 && cond2;
        };

        return ConvexSet<T>(ell, tst, intersect);
    }

    template <typename A, typename B, typename C>
    A first(std::tuple<A, B, C> t)
    {
        return std::get<0>(t);
    }
}
