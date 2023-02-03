#include "stepComp.h"
#include "types.h"

template <typename T>
std::ostream &operator<<(std::ostream &os, const StepComp<T> &s)
{
    if (s.is_done())
    {
        os << "Done(" << s.value() << ")";
    }
    else
    {
        os << "Incomplete(speed = " << s.speed() << ")";
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
Maybe<T> StepComp<T>::get_result() const
{
    return (this->done_) ? this->value_ : Maybe<T>();
}

template <typename T>
StepComp<StepComp<T>> StepComp<T>::subtask(int n) const
{
    if (n <= 0 || this->done_)
    {
        return StepComp<StepComp<T>>(*this);
    }
    StepComp<T> copy = *this;
    return StepComp<StepComp<T>>([copy, n]()
                                 { return copy.untick().subtask(n - 1); });
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
Maybe<T> StepComp<T>::run_bounded(int n) const
{
    return this->forward(n).get_result();
}

namespace stepcomp
{
    /**
     * @brief Mimics the >>= operator from Haskell.
     */
    template <typename T, typename U>
    StepComp<U> bind(StepComp<T> sc, std::function<StepComp<U>(T)> g)
    {
        if (sc.is_done())
        {
            return g(sc.value());
        }
        return StepComp<U>([sc, g]()
                           { return bind(sc.untick(), g); });
    }

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
            current = StepComp<T>([current]()
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
        return StepComp<Either<T1, T2>>([c1, c2]()
                                        { return parallel<A, B>(c1.untick(), c2.untick()); });
    }

    template <typename T>
    StepComp<T> diverge()
    {
        return StepComp<T>([]()
                           { return diverge<T>(); });
    }

    template <typename T>
    StepComp<T> parallel_first(StepComp<T> c1, StepComp<T> c2)
    {
        if (c1.is_done())
        {
            return StepComp<T>(c1.value());
        }
        if (c2.is_done())
        {
            return StepComp<T>(c2.value());
        }
        return StepComp<T>([c1, c2]()
                           { return parallel_first(c1.untick(), c2.untick()); });
    }

    template <typename A, typename B>
    StepComp<Maybe<std::tuple<A, B>>> parallel_maybe(StepComp<Maybe<A>> c1, StepComp<Maybe<B>> c2)
    {
        if (c1.is_done() && c2.is_done())
        {
            if (!c1.value().has_value() || !c2.value().has_value())
            {
                return StepComp<Maybe<std::tuple<A, B>>>(std::nullopt);
            }
            return StepComp(Maybe<std::tuple<A, B>>(std::make_tuple(c1.value().value(), c2.value().value())));
        }
        if (c1.is_done())
        {
            Maybe<A> m1 = c1.value();
            if (!m1.has_value())
            {
                return StepComp<Maybe<std::tuple<A, B>>>(std::nullopt);
            }
            return StepComp<Maybe<std::tuple<A, B>>>([c1, c2]()
                                                     { return parallel_maybe(c1, c2.untick()); });
        }
        if (c2.is_done())
        {
            Maybe<B> m2 = c2.value();
            if (!m2.has_value())
            {
                return StepComp<Maybe<std::tuple<A, B>>>(std::nullopt);
            }
            return StepComp<Maybe<std::tuple<A, B>>>([c1, c2]()
                                                     { return parallel_maybe(c1.untick(), c2); });
        }
        return StepComp<Maybe<std::tuple<A, B>>>([c1, c2]()
                                                 { return parallel_maybe(c1.untick(), c2.untick()); });
    }

    template <typename T>
    StepComp<Maybe<List<T>>> parallel_list_maybe(List<StepComp<Maybe<T>>> steps)
    {
        if (steps.size() == 0)
        {
            return Maybe<List<T>>(List<T>{});
        }
        bool all_done = true;
        for (StepComp<Maybe<T>> sc : steps)
        {
            // If we get any nullopt results, return a nullopt StepComp right away.
            if (sc.is_done() && !sc.value().has_value())
            {
                return StepComp<Maybe<List<T>>>(std::nullopt);
            }
            if (!sc.is_done())
            {
                all_done = false;
            }
        }
        if (all_done)
        {
            List<T> results;
            std::transform(steps.begin(), steps.end(), std::back_inserter(results), [](StepComp<Maybe<T>> sc)
                           { return sc.value().value(); });
            return StepComp(Maybe<List<T>>(results));
        }
        List<StepComp<Maybe<T>>> next;
        std::transform(steps.begin(), steps.end(), std::back_inserter(next), [](StepComp<Maybe<T>> sc)
                       { return sc.untick(); });
        return StepComp<Maybe<List<T>>>([next]()
                                        { return parallel_list_maybe(next); });
    }
}