#include "types.h"
#include "ring.h"
#include <gmpxx.h>
#include <random>

namespace utils
{
    int to_int(Integer n)
    {
        return ring::mpzToInt(n);
    }

    unsigned long to_unsigned_long(Integer n)
    {
        assert(n >= 0);
        assert(n.fits_ulong_p());
        return static_cast<unsigned long>(n.get_ui());
    }

    template <typename T>
    List<T> tail(List<T> x)
    {
        if (x.size() == 0)
        {
            throw std::invalid_argument("Can't take the tail of an empty list");
        }
        return List<T>(x.begin() + 1, x.end());
    }

    template <typename T>
    List<T> concat(List<T> x, List<T> y)
    {
        List<T> result = x;
        result.insert(result.end(), y.begin(), y.end());
        return result;
    }

    Integer div(Integer x, Integer y)
    {
        Integer q;
        mpz_fdiv_q(q.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
        return q;
    }

    Integer mod(Integer x, Integer y)
    {
        Integer r;
        mpz_fdiv_r(r.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
        return r;
    }

    /**
     * Uniform random integer. Both endpoints are inclusive.
     */
    int randint(int lo, int hi)
    {
        assert(lo <= hi);
        std::random_device r;
        std::default_random_engine e1(r());
        std::uniform_int_distribution<int> gen(lo, hi);
        return gen(e1);
    }

    gmp_randclass rr(gmp_randinit_default);
    bool seeded = false;

    /**
     * Uniform random integer. Both endpoints are inclusive.
     */
    Integer randint(Integer lo, Integer hi)
    {
        assert(lo <= hi);

        if (!seeded)
        {
            rr.seed((unsigned long)(time(NULL)));
            seeded = true;
        }

        Integer n = hi - lo + 1;
        Integer rand = rr.get_z_range(n);
        return lo + rand;
    }
}