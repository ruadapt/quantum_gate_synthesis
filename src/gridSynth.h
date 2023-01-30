#pragma once
#include "ring.h"
#include "types.h"
#include "gridproblems.h"


namespace gridsynth
{
    enum DStatus
    {
        Success,
        Fail,
        Timeout,
    };
    
    template <typename T>
    ConvexSet<T> epsilon_region(T epsilon, T theta);

    template <typename T>
    ConvexSet<T> epsilon_region_scaled(DRootTwo s, T epsilon, T theta);

    template <typename T>
    std::tuple<U2<DOmega>, Maybe<double>, List<std::tuple<DOmega, Integer, DStatus>>> gridsynth_internal(
        T prec, T theta, int effort);

    template <typename A, typename B, typename C>
    A first(std::tuple<A, B, C> t);
}

#include "gridSynth.cpp"