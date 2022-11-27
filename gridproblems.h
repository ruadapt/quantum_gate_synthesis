#pragma once
#include "types.h"
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
    Ellipse(Operator<T> op, Point<T> p) : op_{op}, p_{p}
    {
    }
    Operator<T> op() const;
    Point<T> p() const;

private:
    Operator<T> op_;
    Point<T> p_;
};

using CharFun = std::function<bool(Point<DRootTwo>)>;

template <typename T>
using LineIntersector = std::function<std::optional<std::tuple<T, T>>(Point<DRootTwo>, Point<DRootTwo>)>;

template <typename T>
class ConvexSet
{
public:
    ConvexSet(Ellipse<T> el, CharFun test, LineIntersector<T> intersect) : el_{el}, test_{test}, intersect_{intersect}
    {
    }
    Ellipse<T> el() const;
    CharFun test() const;
    LineIntersector<T> intersect() const;

private:
    Ellipse<T> el_;
    CharFun test_;
    LineIntersector<T> intersect_;
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
    T iprod(Point<T> p1, Point<T> p2);

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

    template <typename T>
    ConvexSet<T> unitDisk();
}

#include "gridproblems.cpp"