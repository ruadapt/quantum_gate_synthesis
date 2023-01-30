#pragma once
#include "types.h"
#include "ring.h"

template <typename T>
T to_real(Rational r);

template <typename T>
T to_real(Integer n);

template <typename T>
T to_real(int n);

template <typename T>
T to_real(Real d);

#include "toReal.cpp"