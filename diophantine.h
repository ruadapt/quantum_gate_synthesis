#pragma once
#include "ring.h"
#include "stepComp.h"

namespace diophantine
{
    StepComp<Maybe<ZOmega>> diophantine(ZRootTwo xi);

    StepComp<Maybe<DOmega>> diophantine_dyadic(DRootTwo xi);

    StepComp<Maybe<ZRootTwo>> diophantine_associate(ZRootTwo xi);

    StepComp<Integer> find_factor(Integer n);

    template <typename T>
    List<std::tuple<T, List<std::tuple<T, Integer>>>> relatively_prime_factors(T a, T b);

    Integer power_mod(Integer a, Integer k, Integer n);

    StepComp<Integer> root_of_negative_one(Integer n);

    StepComp<Integer> root_mod(Integer n, Integer a);
}

#include "diophantine.cpp"