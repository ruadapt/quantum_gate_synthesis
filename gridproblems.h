#pragma once
#include "ring.h"
#include <tuple>

namespace gridprob
{
    template <typename T>
    T lambda();

    template <typename T>
    T lambdaInv();

    template <typename T>
    bool within(T x, T low, T high);

    template <typename T>
    std::tuple<Integer, T> floorlog(T b, T x);

    template <typename T>
    std::vector<ZRootTwo> gridpointsInternal(T x0, T x1, T y0, T y1);

    template <typename T>
    std::vector<ZRootTwo> gridpoints(T x0, T x1, T y0, T y1);

    template <typename T>
    std::vector<DRootTwo> gridpointsScaled(T x0, T x1, T y0, T y1, Integer k);
}

#include "gridproblems.cpp"