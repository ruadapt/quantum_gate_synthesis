#pragma once
#include "ring.h"
#include "stepComp.h"

namespace diophantine
{
    StepComp<Maybe<ZOmega>> diophantine(ZRootTwo xi);

    StepComp<Maybe<DOmega>> diophantine_dyadic(DRootTwo xi);

    StepComp<Maybe<ZRootTwo>> diophantine_associate(ZRootTwo xi);

    StepComp<Integer> find_factor(Integer n);

    std::tuple<Integer, List<std::tuple<Integer, Integer>>> relatively_prime_factors(Integer a, Integer b);

    Integer power_mod(Integer a, Integer k, Integer n);

    StepComp<Integer> root_of_negative_one(Integer n);

    StepComp<Integer> root_mod(Integer n, Integer a);
}

#include "diophantine.cpp"