#include "ring.h"
#include <vector>

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
}