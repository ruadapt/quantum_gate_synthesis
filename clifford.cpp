#include "clifford.h"
#include "utils.h"

namespace clifford
{
    Clifford clifford_X()
    {
        return Clifford(0, 1, 0, 0);
    }

    Clifford clifford_Y()
    {
        return Clifford(0, 1, 2, 2);
    }

    Clifford clifford_Z()
    {
        return Clifford(0, 0, 2, 0);
    }

    Clifford clifford_H()
    {
        return Clifford(1, 0, 1, 5);
    }

    Clifford clifford_S()
    {
        return Clifford(0, 0, 1, 0);
    }

    Clifford clifford_SH()
    {
        return clifford_mult(clifford_S(), clifford_H());
    }

    Clifford clifford_E()
    {
        return Clifford(1, 0, 0, 0);
    }

    Clifford clifford_W()
    {
        return Clifford(0, 0, 0, 1);
    }

    template <typename T>
    Clifford to_clifford(T arg)
    {
        static_assert(is_list_v<T>);
        if (arg.empty())
        {
            return clifford_id();
        }
        // TODO convert this to an iterative approach.
        return clifford_mult(to_clifford(arg.front()), to_clifford(utils::tail(arg)));
    }

    template <>
    Clifford to_clifford(Clifford c)
    {
        return c;
    }

    template <>
    Clifford to_clifford(char c)
    {
        switch (c)
        {
        case 'E':
            return clifford_E();
        case 'X':
            return clifford_X();
        case 'S':
            return clifford_S();
        case 'W':
            return clifford_W();
        case 'I':
            return clifford_id();
        case 'i':
            return Clifford(0, 0, 0, 2);
        case '-':
            return Clifford(0, 0, 0, 4);
        case 'H':
            return clifford_H();
        case 'Y':
            return clifford_Y();
        case 'Z':
            return clifford_Z();
        default:
            throw std::invalid_argument("Invalid character for to_clifford");
        }
    }

    /**
     * Unlike in Haskell, strings in C++ aren't just a list of chars, so we have to
     * include this for strings to be convertible to Clifford.
     */
    template <>
    Clifford to_clifford(std::string s)
    {
        List<char> lst(s.begin(), s.end());
        return to_clifford(lst);
    }

    template <>
    Clifford to_clifford(Axis ax)
    {
        switch (ax)
        {
        case Axis_I:
            return to_clifford('I');
        case Axis_H:
            return to_clifford('H');
        case Axis_SH:
            return to_clifford(std::string("SH"));
        }
    }

    template <typename T>
    Tup4<int> clifford_decompose(T m)
    {
        Clifford c = to_clifford<T>(m);
        return {c.a(), c.b(), c.c(), c.d()};
    }

    template <typename T>
    std::tuple<Axis, int, int, int> clifford_decompose_coset(T u)
    {
        Clifford op = to_clifford(u);
        if (op.a() == 0)
        {
            return {Axis_I, op.b(), op.c(), op.d()};
        }
        if (op.a() == 1)
        {
            Clifford op2 = clifford_mult(clifford_inv('H'), op);
            // TODO if op2.a() != 0, is this an error or should we just move on to the
            // next case.
            assert(op2.a() == 0);
            return {Axis_H, op2.b(), op2.c(), op2.d()};
        }
        if (op.a() == 2)
        {
            Clifford op2 = clifford_mult(clifford_inv(std::string("SH")), op);
            assert(op2.a() == 0); // TODO see comment above about this
            return {Axis_SH, op2.b(), op2.c(), op2.d()};
        }
        throw std::invalid_argument("Invalid input to clifford_decompose_coset");
    }

    Clifford clifford_id()
    {
        return Clifford(0, 0, 0, 0);
    }

    Clifford clifford_mult(Clifford u1, Clifford u2)
    {
        int a3, b3, c3, d3;
        std::tie(a3, b3, c3, d3) = conj3(u1.b(), u1.c(), u2.a());
        int c4, d4;
        std::tie(c4, d4) = conj2(c3, u2.b());
        int a = utils::mod(u1.a() + a3, 3);
        int b = utils::mod(b3 + u2.b(), 2);
        int c = utils::mod(c4 + u2.c(), 4);
        int d = utils::mod(d4 + d3 + u1.d() + u2.d(), 8);
        return Clifford(a, b, c, d);
    }

    template <typename T>
    Clifford clifford_inv(T op)
    {
        Clifford c = to_clifford(op);
        int a2, b2, c2, d2;
        std::tie(a2, b2, c2, d2) = cinv(c.a(), c.b(), c.c());
        int d3 = utils::mod(d2 - c.d(), 8);
        return Clifford(a2, b2, c2, d3);
    }

    std::tuple<Axis, Clifford> clifford_tconj(Clifford u)
    {
        Axis k;
        int c2, d2;
        std::tie(k, c2, d2) = tconj(u.a(), u.b());
        int c = utils::mod(c2 + u.c(), 4);
        int d = utils::mod(d2 + u.d(), 8);
        return {k, Clifford(0, u.b(), c, d)};
    }

    Pair<int> conj2(int a, int b)
    {
        if (a == 0 && b == 0)
            return {0, 0};
        if (a == 0 && b == 1)
            return {0, 0};
        if (a == 1 && b == 0)
            return {1, 0};
        if (a == 1 && b == 1)
            return {3, 2};
        if (a == 2 && b == 0)
            return {2, 0};
        if (a == 2 && b == 1)
            return {2, 4};
        if (a == 3 && b == 0)
            return {3, 0};
        if (a == 3 && b == 1)
            return {1, 6};
        throw std::runtime_error("This code should not be reached");
    }

    Tup4<int> conj3(int a, int b, int c)
    {
        if (a == 0 && b == 0 && c == 0)
            return {0, 0, 0, 0};
        if (a == 0 && b == 0 && c == 1)
            return {1, 0, 0, 0};
        if (a == 0 && b == 0 && c == 2)
            return {2, 0, 0, 0};
        if (a == 0 && b == 1 && c == 0)
            return {0, 0, 1, 0};
        if (a == 0 && b == 1 && c == 1)
            return {2, 0, 3, 6};
        if (a == 0 && b == 1 && c == 2)
            return {1, 1, 3, 4};
        if (a == 0 && b == 2 && c == 0)
            return {0, 0, 2, 0};
        if (a == 0 && b == 2 && c == 1)
            return {1, 1, 2, 2};
        if (a == 0 && b == 2 && c == 2)
            return {2, 1, 0, 0};
        if (a == 0 && b == 3 && c == 0)
            return {0, 0, 3, 0};
        if (a == 0 && b == 3 && c == 1)
            return {2, 1, 3, 6};
        if (a == 0 && b == 3 && c == 2)
            return {1, 0, 1, 2};
        if (a == 1 && b == 0 && c == 0)
            return {0, 1, 0, 0};
        if (a == 1 && b == 0 && c == 1)
            return {1, 0, 2, 0};
        if (a == 1 && b == 0 && c == 2)
            return {2, 1, 2, 2};
        if (a == 1 && b == 1 && c == 0)
            return {0, 1, 1, 0};
        if (a == 1 && b == 1 && c == 1)
            return {2, 1, 1, 0};
        if (a == 1 && b == 1 && c == 2)
            return {1, 1, 1, 0};
        if (a == 1 && b == 2 && c == 0)
            return {0, 1, 2, 0};
        if (a == 1 && b == 2 && c == 1)
            return {1, 1, 0, 6};
        if (a == 1 && b == 2 && c == 2)
            return {2, 0, 2, 6};
        if (a == 1 && b == 3 && c == 0)
            return {0, 1, 3, 0};
        if (a == 1 && b == 3 && c == 1)
            return {2, 0, 1, 4};
        if (a == 1 && b == 3 && c == 2)
            return {1, 0, 3, 2};
        throw std::runtime_error("This code should not be reached");
    }

    Tup4<int> cinv(int a, int b, int c)
    {
        if (a == 0 && b == 0 && c == 0)
            return {0, 0, 0, 0};
        if (a == 0 && b == 0 && c == 1)
            return {0, 0, 3, 0};
        if (a == 0 && b == 0 && c == 2)
            return {0, 0, 2, 0};
        if (a == 0 && b == 0 && c == 3)
            return {0, 0, 1, 0};
        if (a == 0 && b == 1 && c == 0)
            return {0, 1, 0, 0};
        if (a == 0 && b == 1 && c == 1)
            return {0, 1, 1, 6};
        if (a == 0 && b == 1 && c == 2)
            return {0, 1, 2, 4};
        if (a == 0 && b == 1 && c == 3)
            return {0, 1, 3, 2};
        if (a == 1 && b == 0 && c == 0)
            return {2, 0, 0, 0};
        if (a == 1 && b == 0 && c == 1)
            return {1, 0, 1, 2};
        if (a == 1 && b == 0 && c == 2)
            return {2, 1, 0, 0};
        if (a == 1 && b == 0 && c == 3)
            return {1, 1, 3, 4};
        if (a == 1 && b == 1 && c == 0)
            return {2, 1, 2, 2};
        if (a == 1 && b == 1 && c == 1)
            return {1, 1, 1, 6};
        if (a == 1 && b == 1 && c == 2)
            return {2, 0, 2, 2};
        if (a == 1 && b == 1 && c == 3)
            return {1, 0, 3, 4};
        if (a == 2 && b == 0 && c == 0)
            return {1, 0, 0, 0};
        if (a == 2 && b == 0 && c == 1)
            return {2, 1, 3, 6};
        if (a == 2 && b == 0 && c == 2)
            return {1, 1, 2, 2};
        if (a == 2 && b == 0 && c == 3)
            return {2, 0, 3, 6};
        if (a == 2 && b == 1 && c == 0)
            return {1, 0, 2, 0};
        if (a == 2 && b == 1 && c == 1)
            return {2, 1, 1, 6};
        if (a == 2 && b == 1 && c == 2)
            return {1, 1, 0, 2};
        if (a == 2 && b == 1 && c == 3)
            return {2, 0, 1, 6};
        throw std::runtime_error("This code should not be reached");
    }

    std::tuple<Axis, int, int> tconj(int a, int b)
    {
        if (a == 0 && b == 0)
            return {Axis_I, 0, 0};
        if (a == 0 && b == 1)
            return {Axis_I, 1, 7};
        if (a == 1 && b == 0)
            return {Axis_H, 3, 3};
        if (a == 1 && b == 1)
            return {Axis_H, 2, 0};
        if (a == 2 && b == 0)
            return {Axis_SH, 0, 5};
        if (a == 2 && b == 1)
            return {Axis_SH, 1, 4};
        throw std::runtime_error("This code should not be reached");
    }
}