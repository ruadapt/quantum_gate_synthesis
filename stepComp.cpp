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
        os << "StepComp(Incomplete, speed = " << s.speed() << ")";
    }
    return os;
}

template <typename T>
bool StepComp<T>::is_done() const { return done_; }

template <typename T>
int StepComp<T>::speed() const { return speed_; }

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
int StepComp<T>::count() const
{
    return count_;
}

template <typename T>
StepComp<T> StepComp<T>::set_count(int count) const
{
    StepComp<T> result = *this;
    result.count_ = count;
    return result;
}

template <typename T>
void StepComp<T>::reset_count() const
{
    count_ = 0;
}

template <typename T>
StepComp<T> StepComp<T>::untick() const
{
    if (this->is_done())
    {
        return StepComp<T>(this->value(), this->speed_).set_count(count_);
    }
    int c = this->count();
    StepComp<T> current = *this;
    for (int i = 0; i < this->speed_; i++)
    {
        if (current.is_done())
        {
            break;
        }
        current = current.comp_();
        c++;
    }
    // The speed always stays the same after unticking, regardless of the contents of the
    // inner StepComps.
    current.speed_ = speed_;
    return current.set_count(c);
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

template <typename T>
std::optional<T> StepComp<T>::get_result() const
{
    return (this->done_) ? this->value_ : std::nullopt;
}

template <typename T>
StepComp<StepComp<T>> StepComp<T>::subtask(int n) const
{
    if (n <= 0 || this->done_)
    {
        return StepComp<StepComp<T>>(*this);
    }
    // TODO try to make this iterative.
    return StepComp<StepComp<T>>([=]()
                                 { this->untick().subtask(n - 1); });
}

template <typename T>
StepComp<T> StepComp<T>::speedup(int n) const
{
    if (this->done_)
    {
        return StepComp<T>(this->value_, this->speed_ * n).set_count(count_);
    }
    return StepComp<T>(this->comp_, this->speed_ * n).set_count(count_);
}

template <typename T>
T StepComp<T>::run() const
{
    StepComp<T> current = *this;
    while (!current.done_)
    {
        current = current.untick();
    }
    return current.value_;
}

template <typename T>
std::optional<T> StepComp<T>::run_bounded(int n)
{
    return this->forward(n).get_result();
}

namespace stepcomp
{
    /**
     * Given a value, put it into a StepComp that needs to run for n steps to
     * finish and yield the value.
     */
    template <typename T>
    StepComp<T> wrap(T value, int n, int speed)
    {
        StepComp<T> current = StepComp<T>(value, speed);
        for (int i = 0; i < n; i++)
        {
            current = StepComp<T>([=]()
                                  { return current; },
                                  speed);
        }
        return current;
    }

    template <typename A, typename B>
    StepComp<Either<std::tuple<A, StepComp<B>>, std::tuple<StepComp<A>, B>>> parallel(StepComp<A> c1, StepComp<B> c2)
    {
        using T1 = std::tuple<A, StepComp<B>>;
        using T2 = std::tuple<StepComp<A>, B>;
        if (c1.is_done())
        {
            return StepComp(Either<T1, T2>(std::make_tuple(c1.value(), c2)));
        }
        if (c2.is_done())
        {
            return StepComp(Either<T1, T2>(std::make_tuple(c1, c2.value())));
        }
        return StepComp<Either<T1, T2>>([=]()
                                        { return parallel<A, B>(c1.untick(), c2.untick()); });
    }

    template <typename T>
    StepComp<T> diverge()
    {
        return StepComp<T>([]()
                           { return diverge<T>(); });
    }
}