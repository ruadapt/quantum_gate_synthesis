#include "multiQubitSynthesis.h"
#include "ring.h"

namespace multi_qubit_synthesis
{
    Z2 residue(Integer n)
    {
        return ring::parity(n);
    }

    template <typename A, typename B>
    Omega<B> residue(Omega<A> o)
    {
        return Omega<B>(residue(o.A()), reside(o.b()), reside(o.c()), reside(o.d()));
    }

    template <typename A, typename B>
    Omega<B> residue(RootTwo<A> r)
    {
        return RootTwo<B>(residue(r.A()), reside(r.b()));
    }
}