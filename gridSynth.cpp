#include "gridSynth.h"
#include "quadratic.h"
#include "ring.h"
#include "diophantine.h"

namespace gridsynth
{
    namespace dio = diophantine;
    namespace gp = gridprob;
    namespace sc = stepcomp;
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

        auto intersect = [d, z](Point<DRootTwo> p, Point<DRootTwo> v) -> Maybe<Pair<T>>
        {
            DRootTwo a = gp::iprod(v, v);
            DRootTwo b = DRootTwo(2) * gp::iprod(v, p);
            DRootTwo c = gp::iprod(p, p) - 1;
            // std::cout << "a = " << a << std::endl;
            // std::cout << "b = " << b << std::endl;
            // std::cout << "c = " << c << std::endl;
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
            // std::cout << "t0 = " << t0 << std::endl;
            // std::cout << "t1 = " << t1 << std::endl;
            // std::cout << "t2 = " << t2 << std::endl;
            if (vz > 0)
            {
                return Maybe<Pair<T>>(Pair<T>{std::max(t0, t2), t1});
            }
            return Maybe<Pair<T>>(Pair<T>{t0, std::min(t1, t2)});
        };

        auto tst = [d, zx, zy](Point<DRootTwo> p) -> bool
        {
            DRootTwo x, y;
            std::tie(x, y) = p;
            bool cond1 = x * x + y * y <= DRootTwo(1);
            bool cond2 = zx * ring::fromDRootTwo<T>(x) + zy * ring::fromDRootTwo<T>(y) >= d;
            return cond1 && cond2;
        };

        return ConvexSet<T>(ell, tst, intersect);
    }

    template <typename T>
    std::tuple<U2<DOmega>, Maybe<double>, List<std::tuple<DOmega, Integer, DStatus>>> gridsynth_internal(
        T prec, T theta, int effort)
    {
        T epsilon = bmp::pow(2, -prec);
        ConvexSet<T> region = epsilon_region(epsilon, theta);
        std::function<List<DOmega>(Integer)> raw_candidates = gp::gridpoints2_increasing(region, gp::unitDisk<T>());

        auto tcount = [](Integer k) -> Integer
        {
            return (k > 0) ? (2 * k - 2) : 0_mpz;
        };

        List<std::tuple<DOmega, Integer, DStatus>> candidate_info;
        Integer k = 0;
        while (true)
        {
            List<DOmega> candidates = raw_candidates(k);
            for (DOmega u : candidates)
            {
                Integer tc = tcount(k);
                RootTwo<ZDyadic> xi = ring::real(DOmega(1) - u.adj() * u);
                Maybe<Maybe<DOmega>> answer_t = dio::diophantine_dyadic(xi).run_bounded(effort);
                if (!answer_t.has_value())
                {
                    candidate_info.push_back({u, tc, Timeout});
                    continue;
                }
                if (!answer_t.value().has_value())
                {
                    candidate_info.push_back({u, tc, Fail});
                    continue;
                }
                candidate_info.push_back({u, tc, Success});
                DOmega t = answer_t.value().value();
                DOmega omega = ring::omega<DOmega>();
                U2<DOmega> uU;
                if (ring::denomExp(u + t) < ring::denomExp(u + omega * t))
                {
                    uU = matrix2x2(u, -t.adj(), t, u.adj());
                }
                else
                {
                    uU = matrix2x2(u, -((omega * t).adj()), omega * t, u.adj());
                }
                return { uU, 1.0, candidate_info }; 
                // TODO add error calculation
            }
            k++;
        }
    }

    template <typename A, typename B, typename C>
    A first(std::tuple<A, B, C> t)
    {
        return std::get<0>(t);
    }
}
