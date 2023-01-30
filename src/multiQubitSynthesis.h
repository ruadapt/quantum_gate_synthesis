/** \file multiQubitSynthesis.h
 */
#pragma once
#include "types.h"

using Index = size_t;

enum TwoLevelType
{
    TL_X,
    TL_H,
    TL_T,
    TL_omega
};

enum ResidueType
{
    RT_0000,
    RT_0001,
    RT_1010
};

class TwoLevel
{
public:
    /**
     * Do not call this constructor directly. Instead, use make_TL_X, make_TL_H,
     * make_TL_T, or make_TL_omega.
     */
    TwoLevel(TwoLevelType type, int pow, Index i1, Index i2) : type_{type}, pow_{pow}, i1_{i1}, i2_{i2}
    {
    }
    bool operator==(const TwoLevel &tl) const
    {
        return (type_ == tl.type_) && (pow_ == tl.pow_) && (i1_ == tl.i1_) && (i2_ == tl.i2_);
    }
    TwoLevelType type() const { return type_; }
    int pow() const
    {
        if (type_ != TL_T && type_ != TL_omega)
        {
            throw std::runtime_error("pow is only defined for TL_T and TL_omega");
        }
        return pow_;
    }
    Index i1() const { return i1_; }
    Index i2() const
    {
        if (type_ == TL_omega)
        {
            throw std::runtime_error("i2 is not defined for TL_omega");
        }
        return i2_;
    }

private:
    TwoLevelType type_;
    int pow_;
    Index i1_;
    Index i2_;
};

std::ostream &operator<<(std::ostream &os, const TwoLevel &tl)
{
    switch (tl.type())
    {
            case TL_X:
            {
                os << "TL_X(" << tl.i1() << ", " << tl.i2() << ")";
                break;
            }
            case TL_H:
            {
                os << "TL_H(" << tl.i1() << ", " << tl.i2() << ")";
                break;
            }
            case TL_T:
            {
                os << "TL_T(" << tl.pow() << ", " << tl.i1() << ", " << tl.i2() << ")";
                break;
            }
            case TL_omega:
            {
                os << "TL_omega(" << tl.pow() << ", " << tl.i1() << ")";
                break;
            }
    }
    return os;
}


TwoLevel make_TL_X(Index i1, Index i2)
{
    return TwoLevel(TL_X, 0, i1, i2);
}

TwoLevel make_TL_H(Index i1, Index i2)
{
    return TwoLevel(TL_H, 0, i1, i2);
}

TwoLevel make_TL_T(int pow, Index i1, Index i2)
{
    return TwoLevel(TL_T, pow, i1, i2);
}

TwoLevel make_TL_omega(int pow, Index i1)
{
    return TwoLevel(TL_omega, pow, i1, 0);
}

#include "multiQubitSynthesis.cpp"