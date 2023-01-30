#pragma once

namespace euclidean_domain
{
    Integer rounddiv(Integer x, Integer y);

    std::tuple<Integer, Integer> divMod(Integer x, Integer y);
}

#include "euclideanDomain.cpp"