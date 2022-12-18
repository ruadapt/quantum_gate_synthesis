#include "stepComp.h"

template <typename T>
StepComp<T>::StepComp(T value)
{
    value_ = value;
    done_ = true;
}

template <typename T>
StepComp<T>::StepComp(std::function<StepComp<T>()> comp)
{
    comp_ = comp;
    done_ = false;
}

template <typename T>
bool StepComp<T>::done() const { return done_; }

template <typename T>
T StepComp<T>::value() const { return value_; }

template <typename T>
std::function<StepComp<T>()> StepComp<T>::comp() const { return comp_; }

template <typename T>
std::ostream &operator<<(std::ostream &os, const StepComp<T> &s)
{
    if (s.done())
    {
        os << "StepComp(" << s.value() << ")";
    }
    else
    {
        os << "StepComp(Incomplete)";
    }
    return os;
}