#include "cliffordT.h"
#include "clifford.h"
#include "multiQubitSynthesis.h"
#include <boost/algorithm/string/join.hpp>

namespace clifford_t
{
    namespace cliff = clifford;
    namespace mqs = multi_qubit_synthesis;

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
    List<Gate> to_gates(T other_form)
    {
        static_assert(is_list_v<T>, "This implementation only works for lists");
        using I = typename T::value_type;
        List<List<Gate>> all_gates = utils::map<I, List<Gate>>(to_gates<I>, other_form);
        return utils::concat(all_gates);
    }

    template <>
    List<Gate> to_gates(Gate x)
    {
        return List<Gate>{x};
    }

    template <>
    List<Gate> to_gates(char ch);

    template <>
    List<Gate> to_gates(Axis ax)
    {
        switch(ax)
        {
        case Axis_I:
            return List<Gate>{};
        case Axis_H:
            return List<Gate>{H};
        case Axis_SH:
            return List<Gate>{S, H};
        }
    }

    template <>
    List<Gate> to_gates(Clifford op)
    {
        Axis k;
        int b, c, d;
        std::tie(k, b, c, d) = cliff::clifford_decompose_coset(op);
        List<Gate> as = to_gates(k);
        List<Gate> xs(size_t(b), X);
        List<Gate> ss(size_t(c), S);
        List<Gate> ws(size_t(d), W);
        return utils::concat(List<List<Gate>>{as, xs, ss, ws});
    }

    // TODO test
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

    template <>
    List<Gate> to_gates(Syllables s)
    {
        List<Gate> reverse_result{};
        for (TrailingSyllable t : s.tail())
        {
            switch (t)
            {
            case HT:
            {
                reverse_result.push_back(T);
                reverse_result.push_back(H);
                break;
            }
            case SHT:
            {
                reverse_result.push_back(T);
                reverse_result.push_back(H);
                reverse_result.push_back(S);
                break;
            }
            }
        }
        if (s.lead() == S_T)
        {
            reverse_result.push_back(T);
        }
        // Reverse before returning.
        std::reverse(reverse_result.begin(), reverse_result.end());
        return reverse_result;
    }

    template <>
    List<Gate> to_gates(NormalForm n)
    {
        return utils::concat(to_gates(n.ts()), to_gates(n.c()));
    }

    template <typename T>
    T from_gates(List<Gate> gates);

    template <>
    std::string from_gates(List<Gate> gates);

    template <>
    List<Gate> from_gates(List<Gate> gates)
    {
        return gates;
    }

    List<Gate> invert_gates(List<Gate> gs);

    template <typename A, typename B>
    B convert(A arg)
    {
        return from_gates<B>(to_gates<A>(arg));
    }

    NormalForm normalform_append(NormalForm n, Gate G)
    {
        switch (G)
        {
        case X:
            return NormalForm(n.ts(), cliff::clifford_mult(n.c(), cliff::clifford_X()));
        case Y:
            return NormalForm(n.ts(), cliff::clifford_mult(n.c(), cliff::clifford_Y()));
        case Z:
            return NormalForm(n.ts(), cliff::clifford_mult(n.c(), cliff::clifford_Z()));
        case H:
            return NormalForm(n.ts(), cliff::clifford_mult(n.c(), cliff::clifford_H()));
        case S:
            return NormalForm(n.ts(), cliff::clifford_mult(n.c(), cliff::clifford_S()));
        case E:
            return NormalForm(n.ts(), cliff::clifford_mult(n.c(), cliff::clifford_E()));
        case W:
            return NormalForm(n.ts(), cliff::clifford_mult(n.c(), cliff::clifford_W()));
        case T:
        {
            Axis k;
            Clifford c_prime;
            std::tie(k, c_prime) = cliff::clifford_tconj(n.c());
            switch (k)
            {
            case Axis_H:
                return NormalForm(cons_SApp_HT(n.ts()), c_prime);
            case Axis_SH:
                return NormalForm(cons_SApp_SHT(n.ts()), c_prime);
            default:
            {
                switch (n.ts().cons())
                {
                case CONS_S_I:
                    return NormalForm(cons_S_T(), c_prime);
                case CONS_S_T:
                    return NormalForm(cons_S_I(), cliff::clifford_mult(cliff::clifford_S(), c_prime));
                case CONS_SApp_HT:
                {
                    Clifford clifford_HS = cliff::to_clifford(std::string("HS"));
                    Syllables ts_prime = Syllables(n.ts().lead(), utils::tail(n.ts().tail()));
                    return NormalForm(ts_prime, cliff::clifford_mult(clifford_HS, c_prime));
                }
                case CONS_SApp_SHT:
                {
                    Clifford clifford_SHS = cliff::to_clifford(std::string("SHS"));
                    Syllables ts_prime = Syllables(n.ts().lead(), utils::tail(n.ts().tail()));
                    return NormalForm(ts_prime, cliff::clifford_mult(clifford_SHS, c_prime));
                }
                }
            }
            }
        }
        }
    }

    NormalForm nf_id()
    {
        return NormalForm(cons_S_I(), cliff::clifford_id());
    }

    template <typename T>
    NormalForm nf_mult(NormalForm a, T b)
    {
        return utils::foldl<Gate, NormalForm>(normalform_append, a, to_gates(b));
    }

    template <typename T>
    NormalForm nf_inv(T gate_convertible)
    {
        return from_gates(invert_gates(to_gates(gate_convertible)));
    }

    template <typename T>
    NormalForm normalize(T gate_convertible)
    {
        return nf_mult(nf_id(), gate_convertible);
    }

    List<Gate> synthesis_u2(U2<DOmega> u)
    {
        return to_gates(normalize(mqs::synthesis_nqubit(u)));
    }
}