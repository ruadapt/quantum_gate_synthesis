#pragma once
#include "types.h"
#include <optional>
#include <tuple>

template <typename T>
class StepComp
{
public:
    StepComp(T value);
    StepComp(std::function<StepComp<T>()> comp);
    bool is_done() const;
    T value() const;
    std::function<StepComp<T>()> comp() const;
    StepComp<T> untick() const;
    StepComp<T> forward(int n) const;
    std::optional<T> get_result() const;
    StepComp<StepComp<T>> subtask(int n) const;
    StepComp<T> speedup(int n) const;
    StepComp<std::tuple<T, int>> with_counter() const;
    T run() const;
    std::tuple<T, int> run_with_steps() const;
    std::optional<T> run_bounded(int n);

private:
    bool done_;
    T value_;
    std::function<StepComp<T>()> comp_;
};

namespace stepcomp
{
    template <typename A, typename B>
    Either<std::tuple<A, StepComp<B>>, std::tuple<StepComp<A>, B>> parallel(StepComp<A> c1, StepComp<B> c2);

    template <typename T>
    StepComp<T> diverge();

    template <typename T>
    StepComp<T> parallel_first(StepComp<T> c1, StepComp<T> c2);

    template <typename A, typename B>
    StepComp<Maybe<std::tuple<A, B>>> parallel_maybe(StepComp<Maybe<A>> c1, StepComp<Maybe<B>> c2);

    template <typename T>
    StepComp<Maybe<List<T>>> parallel_list_maybe(List<StepComp<Maybe<T>>>);

    template <typename T>
    StepComp<T> wrap(T value, int n);
}

#include "stepComp.cpp"