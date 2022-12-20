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
}