#include "types.h"
#include "ring.h"
#include <gmpxx.h>
#include <random>

namespace utils
{
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
        std::random_device r;
        std::default_random_engine e1(r());
        std::uniform_int_distribution<int> gen(lo, hi);
        return gen(e1);
    }

    /**
     * Uniform random integer. Both endpoints are inclusive.
     */
    Integer randint(Integer lo, Integer hi)
    {
        return Integer(randint(ring::mpzToInt(lo), ring::mpzToInt(hi)));
    }
}