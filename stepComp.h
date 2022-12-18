#pragma once

#include "types.h"

template <typename T>
class StepComp
{
public:
    StepComp(T value);
    StepComp(std::function<StepComp<T>()> comp);
    bool done() const;
    T value() const;
    std::function<StepComp<T>()> comp() const;

private:
    bool done_;
    T value_;
    std::function<StepComp<T>()> comp_;
};

#include "stepComp.cpp"