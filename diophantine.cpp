#include "diophantine.h"
#include "types.h"
#include "utils.h"
#include "stepComp.h"
#include "ring.h"
#include "euclideanDomain.h"

namespace sc = stepcomp;
namespace ed = euclidean_domain;

namespace diophantine
{
    StepComp<Integer> find_factor(Integer n)
    {
        if (ring::even(n) && n > 2)
        {
            return StepComp<Integer>(Integer(2));
        }
        Integer a = utils::randint(1, n - 1);
        auto f = [=](Integer x) -> Integer
        { return utils::mod(x * x + a, n); };
        auto aux = [=](Integer x, Integer y) -> StepComp<Integer>
        {
            auto aux_impl = [=](Integer x, Integer y, auto &aux_ref) -> StepComp<Integer>
            {
                Integer d = ed::euclid_gcd(Integer(x - y), n); // TODO which gcd function to use?
                if (d == 1)
                {
                    return StepComp<Integer>([=]()
                                             { return StepComp<Integer>(aux_ref(f(x), f(f(y)), aux_ref)); });
                }
                if (d == n)
                {
                    return find_factor(n);
                }
                return StepComp<Integer>(d);
            };
            return aux_impl(x, y, aux_impl);
        };
        return StepComp<Integer>([=]()
                                 { return aux(Integer(2), f(Integer(2))); });
    }

    std::tuple<Integer, List<Integer>, List<Pair<Integer>>> relatively_prime_aux2(Integer h, List<Pair<Integer>> flist)
    {
        if (flist.size() == 0)
        {
            return std::make_tuple(1, List<Integer>{}, List<Pair<Integer>>{std::make_tuple(h, 1_mpz)});
        }
        Integer f, k;
        std::tie(f, k) = flist.at(0);
        List<Pair<Integer>> fs = utils::tail(flist);
        if (ed::euclid_associates(h, f))
        {
            Integer u2 = ed::euclid_div(h, f);
            List<Pair<Integer>> new_fs = fs;
            new_fs.insert(new_fs.begin(), std::make_tuple(f, k + 1));
            return std::make_tuple(u2, List<Integer>{}, new_fs);
        }
        Integer d = ed::euclid_gcd(h, f);
        if (ed::is_unit(d))
        {
            Integer u;
            List<Integer> hs;
            List<Pair<Integer>> fs2;
            std::tie(u, hs, fs2) = relatively_prime_aux2(h, fs);
            List<Pair<Integer>> new_fs = fs2;
            new_fs.insert(new_fs.begin(), std::make_tuple(f, k));
            return std::make_tuple(u, hs, new_fs);
        }
        List<Integer> temp1{ed::euclid_div(h, d), d};
        List<Integer> temp2(utils::to_unsigned_long(k), ed::euclid_div(f, d));
        List<Integer> temp3(utils::to_unsigned_long(k), d);
        return std::make_tuple(1_mpz, utils::concat(temp1, utils::concat(temp2, temp3)), fs);
    }

    std::tuple<Integer, List<Pair<Integer>>> relatively_prime_aux(Integer u, List<Integer> vals, List<Pair<Integer>> fs)
    {
        if (vals.size() == 0)
        {
            return std::make_tuple(u, fs);
        }
        Integer h = vals.at(0);
        List<Integer> t = utils::tail(vals);
        if (ed::is_unit(h))
        {
            return relatively_prime_aux(h * u, utils::tail(vals), fs);
        }
        Integer u2;
        List<Integer> hs;
        List<Pair<Integer>> fs2;
        std::tie(u2, hs, fs2) = relatively_prime_aux2(h, fs);
        return relatively_prime_aux(u2 * u, utils::concat(hs, t), fs2);
    }

    std::tuple<Integer, List<std::tuple<Integer, Integer>>> relatively_prime_factors(Integer a, Integer b)
    {
        return relatively_prime_aux(1, List<Integer>{a, b}, List<Pair<Integer>>{});
    }

    Integer power_mod(Integer a, Integer k, Integer n)
    {
        if (k == 0)
        {
            return 1;
        }
        if (k == 1)
        {
            return utils::mod(a, n);
        }
        Integer b = power_mod(a, utils::div(k, 2), n);
        if (ring::even(k))
        {
            return utils::mod(b * b, n);
        }
        return utils::mod(b * b * a, n);
    }

    StepComp<Integer> root_of_negative_one(Integer n)
    {
        return StepComp<Integer>([=]()
                                 {
            Integer b = utils::randint(1, n - 1);
            Integer h = power_mod(b, utils::div(n - 1, 4), n);
            Integer r = utils::mod(h * h, n);
            if (r == n - 1)
            {
                return StepComp<Integer>(h);
            }
            if (r != 1)
            {
                return sc::diverge<Integer>();
            }
            // TODO wrap this in a stepcomp?
            return root_of_negative_one(n); });
    }
}