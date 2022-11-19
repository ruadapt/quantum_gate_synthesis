#pragma once
#include "ring.h"
#include "matrix.h"
#include "gridproblems.h"
#include <tuple>

template <typename T>
using Point = std::tuple<T, T>;

template <typename T>
T fst(Point<T> p);

template <typename T>
T snd(Point<T> p);

template <typename T>
using Operator = Matrix<T, 2, 2>;

template <typename T>
class Ellipse
{
public:
    Ellipse(Operator<T> op, Point<T> p);
    Operator<T> op() const;
    Point<T> p() const;

private:
    Operator<T> op_;
    Point<T> p_;
};

namespace gridprob
{
    template <typename T>
    T lambda();

    template <typename T>
    T lambdaInv();

    template <typename T>
    bool within(T x, T low, T high);

    template <typename T>
    std::tuple<Integer, T> floorlog(T b, T x);

    template <typename T>
    std::vector<ZRootTwo> gridpointsInternal(T x0, T x1, T y0, T y1);

    template <typename T>
    std::vector<ZRootTwo> gridpoints(T x0, T x1, T y0, T y1);

    template <typename T>
    std::vector<DRootTwo> gridpointsScaled(T x0, T x1, T y0, T y1, Integer k);

    template <typename T>
    std::vector<DRootTwo> gridpointsScaledParity(DRootTwo beta, T x0, T x1, T y0, T y1, Integer k);

    template <typename T>
    Point<T> pointFromDRootTwo(Point<DRootTwo>);

    template <typename T>
    Operator<T> makeOperator(T x0, T x1, T x2, T x3);
}

#include "gridproblems.cpp"