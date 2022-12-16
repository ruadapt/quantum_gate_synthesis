#include "euclideanDomain.h"
#include "ring.h"

namespace euclidean_domain
{
    Integer rounddiv(Integer x, Integer y)
    {
        Integer temp1, temp2, result;
        mpz_tdiv_q(temp1.get_mpz_t(), y.get_mpz_t(), Integer(2).get_mpz_t());
        temp2 = x + temp1;
        mpz_fdiv_q(result.get_mpz_t(), temp2.get_mpz_t(), y.get_mpz_t());
        return result;
    }
}