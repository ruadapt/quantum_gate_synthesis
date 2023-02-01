/** \file quadratic.h
 */
#pragma once
#include "types.h"
#include <optional>

template <typename T>
Maybe<std::tuple<Real, Real>> quadratic(T a, T b, T c);

#include "quadratic.cpp"