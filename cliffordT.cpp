#include "cliffordT.h"
#include "multiQubitSynthesis.h"
#include <boost/algorithm/string/join.hpp>

namespace clifford_t
{
    std::string to_string(Gate g)
    {
        switch (g)
        {
        case X:
            return "X";
        case Y:
            return "Y";
        case Z:
            return "Z";
        case H:
            return "H";
        case S:
            return "S";
        case T:
            return "T";
        case E:
            return "E";
        case W:
            return "W";
        }
    }

    std::string to_string(List<Gate> gs)
    {
        auto gate_to_string = [](Gate g) -> std::string
        { return to_string(g); };
        return boost::algorithm::join(utils::map<Gate, std::string>(gate_to_string, gs), "");
    }

    template <typename T>
    List<Gate> to_gates(T other_form);

    template <>
    List<Gate> to_gates(TwoLevel tl)
    {
        if (tl == make_TL_X(0, 1))
        {
            return List<Gate>{X};
        }
        if (tl == make_TL_X(1, 0))
        {
            return List<Gate>{X};
        }
        if (tl == make_TL_H(0, 1))
        {
            return List<Gate>{H};
        }
        if (tl == make_TL_H(1, 0))
        {
            return List<Gate>{X, H, X};
        }
        if (tl.type() == TL_T && tl.i1() == 0 && tl.i2() == 1)
        {
            int k = tl.pow();
            if (utils::mod(k, 2) == 1)
            {
                return utils::concat(List<Gate>{T}, to_gates<TwoLevel>(make_TL_T(k - 1, 0, 1)));
            }
            if (utils::mod(k, 4) == 2)
            {
                return utils::concat(List<Gate>{S}, to_gates<TwoLevel>(make_TL_T(k - 2, 0, 1)));
            }
            if (utils::mod(k, 8) == 4)
            {
                return List<Gate>{Z};
            }
            return List<Gate>{};
        }
        if (tl.type() == TL_T && tl.i1() == 1 && tl.i2() == 0)
        {
            int k = tl.pow();
            return utils::concat(utils::concat(List<Gate>{X}, to_gates<TwoLevel>(make_TL_T(k, 0, 1))), List<Gate>{X});
        }
        if (tl.type() == TL_omega && tl.i1() == 1)
        {
            int k = tl.pow();
            return to_gates<TwoLevel>(make_TL_T(k, 0, 1));
        }
        if (tl.type() == TL_omega && tl.i1() == 0)
        {
            int k = tl.pow();
            return to_gates<TwoLevel>(make_TL_T(k, 1, 0));
        }
        throw std::invalid_argument("Invalid TwoLevel input to to_gates");
    }

    template <typename T>
    T from_gates(List<Gate> gates);

    template <typename A, typename B>
    B convert(A arg)
    {
        return from_gates<B>(to_gates<A>(arg));
    }
}