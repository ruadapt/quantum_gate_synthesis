/** \file diophantine.h
 */
#pragma once
#include "ring.h"
#include "stepComp.h"

namespace diophantine
{
    StepComp<Maybe<ZOmega>> diophantine(ZRootTwo xi);

    StepComp<Maybe<DOmega>> diophantine_dyadic(DRootTwo xi);

    StepComp<Maybe<ZOmega>> diophantine_associate(ZRootTwo xi);

    StepComp<Integer> find_factor(Integer n);

    template <typename T>
    std::tuple<T, List<std::tuple<T, Integer>>> relatively_prime_factors(T a, T b);

    Integer power_mod(Integer a, Integer k, Integer n);

    StepComp<Integer> root_of_negative_one(Integer n);

    StepComp<Integer> root_mod(Integer n, Integer a);

    StepComp<Maybe<ZOmega>> dioph_int_assoc_prime(Integer n);

    StepComp<Maybe<ZOmega>> dioph_int_assoc(Integer n);
    
    StepComp<Maybe<ZOmega>> dioph_int_assoc_powers(List<Pair<Integer>> facs);

    StepComp<Maybe<ZOmega>> dioph_int_assoc_power(Pair<Integer> p);

    StepComp<Maybe<ZOmega>> dioph_zroottwo_selfassociate(ZRootTwo xi);

    StepComp<Maybe<ZOmega>> dioph_zroottwo_assoc_prime(ZRootTwo xi);

    StepComp<Maybe<ZOmega>> dioph_zroottwo_assoc(ZRootTwo xi);

    StepComp<Maybe<ZOmega>> dioph_zroottwo_assoc_powers(List<std::tuple<ZRootTwo, Integer>> facs);

    StepComp<Maybe<ZOmega>> dioph_zroottwo_assoc_power(std::tuple<ZRootTwo, Integer> p);
}

#include "diophantine.cpp"