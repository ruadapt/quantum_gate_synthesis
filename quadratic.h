#pragma once
#include "types.h"
#include <optional>

template <typename T>
Maybe<Pair<Integer>> int_quadratic(T b, T c);

template <typename T, typename U>
Maybe<Pair<U>> quadratic_fixedprec(T a, T b, T c);

template <typename T, typename U>
Maybe<Pair<U>> quadratic(T a, T b, T c);

#include "quadratic.cpp"