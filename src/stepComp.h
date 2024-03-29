/** \file stepComp.h
 */
#pragma once
#include "types.h"
#include <optional>
#include <tuple>

/**
 * To mimic the StepComp class from Haskell, some changes had to be made. A StepComp<T>
 * is created either with a value of T, which means it is a finished StepComp, or with
 * a function that returns a new StepComp, which itself can either be done or return
 * another StepComp, and so on.
 * 
 * The speed_ field represents how many steps to move forward at a time. One step refers
 * to calling the inner computation once and returning the new StepComp (or doing nothing
 * if the StepComp is already done).
 * 
 * The count_ field keeps track of the total number of steps that have happened (it stops
 * counting once the StepComp is finished).
 */
template <typename T>
class StepComp
{
public:
    StepComp() : done_{true}, speed_{1}, value_{T()}
    {
    }
    StepComp(T value, int speed = 1) : done_{true}, speed_{speed}, value_{value}
    {
    }
    StepComp(std::function<StepComp<T>()> comp, int speed = 1) : done_{false}, speed_{speed}, value_{T()}, comp_{comp}
    {
    }
    bool is_done() const;
    int speed() const;
    T value() const;
    std::function<StepComp<T>()> comp() const;
    int count() const;
    StepComp<T> set_count(int count) const;
    void reset_count() const;
    StepComp<T> untick() const;
    StepComp<T> forward(int n) const;
    Maybe<T> get_result() const;
    StepComp<StepComp<T>> subtask(int n) const;
    StepComp<T> speedup(int n) const;
    T run() const;
    Maybe<T> run_bounded(int n) const;

private:
    bool done_;
    int speed_;
    T value_;
    std::function<StepComp<T>()> comp_;
    int count_ = 0;
};

namespace stepcomp
{
    template <typename T, typename U>
    StepComp<U> bind(StepComp<T> sc, std::function<StepComp<U>(T)> g);

    template <typename A, typename B>
    StepComp<Either<std::tuple<A, StepComp<B>>, std::tuple<StepComp<A>, B>>> parallel(StepComp<A> c1, StepComp<B> c2);

    template <typename T>
    StepComp<T> diverge();

    template <typename T>
    StepComp<T> parallel_first(StepComp<T> c1, StepComp<T> c2);

    template <typename A, typename B>
    StepComp<Maybe<std::tuple<A, B>>> parallel_maybe(StepComp<Maybe<A>> c1, StepComp<Maybe<B>> c2);

    template <typename T>
    StepComp<Maybe<List<T>>> parallel_list_maybe(List<StepComp<Maybe<T>>> steps);

    template <typename T>
    StepComp<T> wrap(T value, int n, int speed = 1);
}

#include "stepComp.cpp"