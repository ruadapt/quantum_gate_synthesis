#include "stepComp.h"

template <typename T>
std::ostream &operator<<(std::ostream &os, const StepComp<T> &s)
{
    if (s.is_done())
    {
        os << "StepComp(" << s.value() << ")";
    }
    else
    {
        os << "StepComp(Incomplete)";
    }
    return os;
}

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
bool StepComp<T>::is_done() const { return done_; }

template <typename T>
T StepComp<T>::value() const
{
    if (this->done_)
    {
        return value_;
    }
    throw std::invalid_argument("Cannot retrieve value of a StepComp that isn't done");
}

template <typename T>
std::function<StepComp<T>()> StepComp<T>::comp() const
{
    if (!this->done_)
    {
        return comp_;
    }
    throw std::invalid_argument("StepComp is done");
}

template <typename T>
StepComp<T> StepComp<T>::untick() const
{
    if (this->is_done())
    {
        return StepComp<T>(this->value());
    }
    return this->comp_();
}

template <typename T>
StepComp<T> StepComp<T>::forward(int n) const
{
    if (this->done_)
    {
        return *this;
    }
    StepComp<T> current = *this;
    for (int i = 0; i < n; i++)
    {
        current = current.untick();
    }
    return current;
}

namespace stepcomp
{
    /**
     * Given a value, put it into a StepComp that needs to be unticked n times to
     * finish and yield the value.
     */
    template <typename T>
    StepComp<T> wrap(T value, int n)
    {
        StepComp<T> current = StepComp<T>(value);
        for (int i = 0; i < n; i++)
        {
            current = StepComp<T>([=]()
                                  { return current; });
        }
        return current;
    }
}