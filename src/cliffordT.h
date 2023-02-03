/** \file cliffordT.h
 */
#pragma once
#include "ring.h"
#include "types.h"
#include "utils.h"
#include "clifford.h"
#include "multiQubitSynthesis.h"

enum Gate
{
    X,
    Y,
    Z,
    H,
    S,
    T,
    E,
    W
};

/**
 * @brief The names of the two no-arg Syllables constructors in Haskell.
 */
enum LeadSyllable
{
    S_I,
    S_T
};

/**
 * @brief The names of the Syllables constructors in Haskell that take in an argument.
 */
enum TrailingSyllable
{
    HT,
    SHT
};

/**
 * @brief Represents the four constructors for Syllables in Haskell.
 */
enum SyllableConstructor
{
    CONS_S_I,
    CONS_S_T,
    CONS_SApp_HT,
    CONS_SApp_SHT
};

/**
 * This class is laid out slightly different than the Haskell version. The lead field
 * corresponds to whether the Haskell object would have been originally created with
 * S_I or S_T, and the tail field is the sequence of SApp_HT and SApp_SHT constructors
 * that were applied from left to right. In other words, SApp_HT (SApp_SHT S_I) in Haskell
 * would become a c++ object with lead = S_I and tail = [HT, SHT].
 */
class Syllables
{
public:
    Syllables(LeadSyllable lead, List<TrailingSyllable> tail) : lead_{lead}, tail_{tail}
    {
    }
    LeadSyllable lead() const { return lead_; }
    List<TrailingSyllable> tail() const { return tail_; }
    bool operator==(const Syllables &s) const { return lead_ == s.lead_ && tail_ == s.tail_; }
    /**
     * Yields the constructor that would be used for pattern matching in Haskell.
     */
    SyllableConstructor cons() const
    {
        if (is_S_I())
        {
            return CONS_S_I;
        }
        if (is_S_T())
        {
            return CONS_S_T;
        }
        if (is_SApp_HT())
        {
            return CONS_SApp_HT;
        }
        if (is_SApp_SHT())
        {
            return CONS_SApp_SHT;
        }
        throw std::runtime_error("This code should not be reachable.");
    }

private:
    bool is_S_I() const { return tail_.empty() && lead_ == S_I; };
    bool is_S_T() const { return tail_.empty() && lead_ == S_T; };
    bool is_SApp_HT() const { return !tail_.empty() && tail_.front() == HT; }
    bool is_SApp_SHT() const { return !tail_.empty() && tail_.front() == SHT; }
    LeadSyllable lead_;
    List<TrailingSyllable> tail_;
};

/**
 * @brief Construct an equivalent object to the S_I constructor in Haskell.
 */
Syllables cons_S_I()
{
    return Syllables(S_I, List<TrailingSyllable>{});
}

/**
 * @brief Construct an equivalent object to the S_T constructor in Haskell.
 */
Syllables cons_S_T()
{
    return Syllables(S_T, List<TrailingSyllable>{});
}

/**
 * @brief Construct an equivalent object to the SApp_HT constructor in Haskell.
 */
Syllables cons_SApp_HT(Syllables s)
{
    return Syllables(s.lead(), utils::concat(List<TrailingSyllable>{HT}, s.tail()));
}

/**
 * @brief Construct an equivalent object to the SApp_SHT constructor in Haskell.
 */
Syllables cons_SApp_SHT(Syllables s)
{
    return Syllables(s.lead(), utils::concat(List<TrailingSyllable>{SHT}, s.tail()));
}

class NormalForm
{
public:
    NormalForm(Syllables ts, Clifford c) : ts_{ts}, c_{c}
    {
    }
    bool operator==(const NormalForm &n) const { return ts_ == n.ts_ && c_ == n.c_; }
    Syllables ts() const { return ts_; }
    Clifford c() const { return c_; }

private:
    Syllables ts_;
    Clifford c_;
};

namespace clifford_t
{
    std::string to_string(Gate g);

    std::string to_string(List<Gate> gs);

    template <typename T>
    List<Gate> to_gates(T other_form);

    template <>
    List<Gate> to_gates(Gate x);

    template <>
    List<Gate> to_gates(char ch);

    template <>
    List<Gate> to_gates(Axis ax);

    template <>
    List<Gate> to_gates(Clifford op);

    template <>
    List<Gate> to_gates(TwoLevel tl);

    template <>
    List<Gate> to_gates(NormalForm n);

    template <>
    List<Gate> to_gates(Syllables s);

    template <typename T>
    T from_gates(List<Gate> gates);

    template <>
    std::string from_gates(List<Gate> gates);

    template <>
    List<Gate> from_gates(List<Gate> gates);

    List<Gate> invert_gates(List<Gate> gs);

    template <typename A, typename B>
    B convert(A arg);

    NormalForm normalform_append(NormalForm n, Gate G);

    NormalForm nf_id();

    template <typename T>
    NormalForm nf_mult(NormalForm a, T b);

    template <typename T>
    NormalForm nf_inv(T arg);

    template <typename T>
    NormalForm normalize(T arg);

    List<Gate> synthesis_u2(U2<DOmega> u);
}

#include "cliffordT.cpp"