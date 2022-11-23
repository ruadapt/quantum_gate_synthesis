#pragma once
#include "types.h"
#include "ring.h"

template <typename T>
T toReal(Rational r);

template <typename T>
T toReal(Integer n);

template <typename T>
T toReal(int n);

template <typename T>
T toReal(Real d);

#include "toReal.cpp"