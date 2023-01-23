#pragma once
#include "types.h"
#include "utils.h"
#include "clifford.h"

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

enum LeadSyllable
{
    S_I,
    S_T
};

enum TrailingSyllables
{
    HT,
    SHT
};

/**
 * This class is laid out slightly different than the Haskell version. The lead field
 * corresponds to whether the Haskell object would have been originally created with
 * S_I or S_T, and the tail field is the sequence of SApp_HT and SApp_SHT constructors
 * that were applied from right to left. In other words, SApp_HT (SApp_SHT S_I) in Haskell
 * would become a c++ object with lead = S_I and tail = [HT, SHT].
 */
class Syllables
{
public:
    Syllables(LeadSyllable lead, List<TrailingSyllables> tail) : lead_{lead}, tail_{tail}
    {
    }
    LeadSyllable lead() const { return lead_; }
    List<TrailingSyllables> tail() const { return tail_; }
    bool operator==(const Syllables &s) const { return lead_ == s.lead_ && tail_ == s.tail_; }
    bool is_S_I() const {return tail_.empty() && lead_ == S_I; };
    bool is_S_T() const {return tail_.empty() && lead_ == S_T; };
    bool is_SApp_HT() const {return !tail_.empty() && tail_.front() == HT; }
    bool is_SApp_SHT() const {return !tail_.empty() && tail_.front() == SHT; }

private:
    LeadSyllable lead_;
    List<TrailingSyllables> tail_;
};

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

#include "cliffordT.cpp"