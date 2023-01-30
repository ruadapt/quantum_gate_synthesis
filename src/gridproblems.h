/** \file gridproblems.h
 */
#pragma once
#include "types.h"
#include "ring.h"
#include "matrix.h"
#include "gridproblems.h"
#include <tuple>

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
    CharFun test_fun() const;
    LineIntersector<T> intersect_fun() const;
    bool test(Point<DRootTwo> p) const;
    std::optional<std::tuple<T, T>> intersect(Point<DRootTwo> p1, Point<DRootTwo> p2) const;

private:
    Ellipse<T> el_;
    CharFun test_;
    LineIntersector<T> intersect_;
};

namespace gridprob
{
    double logBase(double b, double x);

    template <typename T>
    T lambda();

    template <typename T>
    T lambda_inv();

    template <typename T>
    bool within(T x, T low, T high);

    template <typename T>
    std::tuple<Integer, T> floorlog(T b, T x);

    template <typename T>
    double logBase_double(T b, T x);

    template <typename T>
    T iprod(Point<T> p1, Point<T> p2);

    template <typename T>
    std::vector<ZRootTwo> gridpoints_internal(T x0, T x1, T y0, T y1);

    template <typename T>
    std::vector<ZRootTwo> gridpoints(T x0, T x1, T y0, T y1);

    template <typename T>
    std::vector<DRootTwo> gridpoints_scaled(T x0, T x1, T y0, T y1, Integer k);

    template <typename T>
    std::vector<DRootTwo> gridpoints_scaled_parity(DRootTwo beta, T x0, T x1, T y0, T y1, Integer k);

    template <typename T>
    Point<T> point_fromDRootTwo(Point<DRootTwo>);

    template <typename T>
    ConvexSet<T> unitdisk();

    template <typename T>
    ConvexSet<T> disk(DRootTwo s);

    template <typename T>
    Operator<T> op_fromDRootTwo(Operator<DRootTwo> op);

    template <typename T>
    std::tuple<T, double> operator_to_bz(Operator<T> op);

    template <typename T>
    T det(Operator<T> op);

    template <typename T>
    T operator_skew(Operator<T> op);

    template <typename T>
    T uprightness(Operator<T> op);

    template <typename T>
    T skew(OperatorPair<T> pair);

    template <typename T>
    double bias(OperatorPair<T> pair);
}

#include "gridproblems.cpp"