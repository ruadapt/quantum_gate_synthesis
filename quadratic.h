#pragma once
#include <optional>

template <typename T>
std::optional<std::tuple<double, double>> quadratic(T a, T b, T c);

#include "quadratic.cpp"