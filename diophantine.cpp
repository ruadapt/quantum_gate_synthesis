#include "diophantine.h"
#include "types.h"
#include "utils.h"
#include "stepComp.h"
#include "ring.h"
#include "gridproblems.h"
#include "euclideanDomain.h"

namespace gp = gridprob;
namespace sc = stepcomp;
namespace ed = euclidean_domain;

namespace diophantine
{
    StepComp<Maybe<ZOmega>> diophantine(ZRootTwo xi)
    {
        if (xi == 0)
        {
            return StepComp(Maybe<ZOmega>(ZOmega(0)));
        }
        if (xi < 0)
        {
            return StepComp(Maybe<ZOmega>());
        }
        if (xi.adj2() < 0)
        {
            return StepComp(Maybe<ZOmega>());
        }
        StepComp<Maybe<ZOmega>> sc = diophantine_associate(xi);
        auto g = [=](Maybe<ZOmega> t_opt) -> StepComp<Maybe<ZOmega>>
        {
            if (!t_opt.has_value())
            {
                return StepComp(Maybe<ZOmega>());
            }
            ZOmega t = t_opt.value();
            ZRootTwo xi_associate = ring::zRootTwoOfZOmega(t.adj() * t);
            ZRootTwo u = ed::euclid_div(xi, xi_associate);
            Maybe<ZRootTwo> root = ring::zRootTwoRoot(u);
            if (!root.has_value())
            {
                return StepComp(Maybe<ZOmega>());
            }
            ZRootTwo v = root.value();
            return StepComp(Maybe<ZOmega>(ring::fromZRootTwo<ZOmega>(v) * t));
        };
        return sc::bind<Maybe<ZOmega>, Maybe<ZOmega>>(sc, g);
    }

    StepComp<Maybe<DOmega>> diophantine_dyadic(DRootTwo xi)
    {
        Integer k = ring::denomExp(xi);
        Integer k2, k3;
        std::tie(k2, k3) = ed::divMod(k, 2_mpz);
        DRootTwo prod = ring::powNonNeg((gp::lambda<DRootTwo>() * ring::rootTwo<DRootTwo>()), k3) * DRootTwo(ring::powNonNeg(2_mpz, k2)) * xi;
        ZRootTwo xi2 = ring::toWhole<DRootTwo, ZRootTwo>(prod);
        StepComp<Maybe<ZOmega>> sc = diophantine(xi2);
        auto g = [=](Maybe<ZOmega> t2_opt) -> StepComp<Maybe<DOmega>>
        {
            if (!t2_opt.has_value())
            {
                return StepComp(Maybe<DOmega>());
            }
            ZOmega t2 = t2_opt.value();
            DOmega u = ring::rootHalf<DOmega>() * (ring::omega<DOmega>() - ring::i<DOmega>());
            DOmega prod = ring::powNonNeg(u, k3) * ring::powNonNeg(ring::rootHalf<DOmega>(), k2) * ring::fromWhole<DOmega, ZOmega>(t2);
            return StepComp(Maybe<DOmega>(prod));
        };
        return sc::bind<Maybe<ZOmega>, Maybe<DOmega>>(sc, g);
    }

    StepComp<Maybe<ZOmega>> diophantine_associate(ZRootTwo xi)
    {
        if (xi == 0)
        {
            return StepComp(Maybe<ZOmega>(ZOmega(0)));
        }
        ZRootTwo d = ed::euclid_gcd(xi, xi.adj2());
        ZRootTwo xi2 = ed::euclid_div(xi, d);
        StepComp<Maybe<std::tuple<ZOmega, ZOmega>>> sc = sc::parallel_maybe(
            dioph_zroottwo_selfassociate(d), dioph_zroottwo_assoc(xi2));
        auto g = [=](Maybe<std::tuple<ZOmega, ZOmega>> res) -> StepComp<Maybe<ZOmega>>
        {
            if (!res.has_value())
            {
                return StepComp(Maybe<ZOmega>());
            }
            ZOmega t1, t2;
            std::tie(t1, t2) = res.value();
            return StepComp(Maybe<ZOmega>(t1 * t2));
        };
        return sc::bind<Maybe<std::tuple<ZOmega, ZOmega>>, Maybe<ZOmega>>(sc, g);
    }

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

    template <typename T>
    std::tuple<T, List<T>, List<std::tuple<T, Integer>>> relatively_prime_aux2(T h, List<std::tuple<T, Integer>> flist)
    {
        if (flist.size() == 0)
        {
            return std::make_tuple(T(1), List<T>{}, List<std::tuple<T, Integer>>{std::make_tuple(h, 1_mpz)});
        }
        T f;
        Integer k;
        std::tie(f, k) = flist.at(0);
        List<std::tuple<T, Integer>> fs = utils::tail(flist);
        if (ed::euclid_associates(h, f))
        {
            T u2 = ed::euclid_div(h, f);
            List<std::tuple<T, Integer>> new_fs = fs;
            new_fs.insert(new_fs.begin(), std::make_tuple(f, k + 1));
            return std::make_tuple(u2, List<T>{}, new_fs);
        }
        T d = ed::euclid_gcd(h, f);
        if (ed::is_unit(d))
        {
            T u;
            List<T> hs;
            List<std::tuple<T, Integer>> fs2;
            std::tie(u, hs, fs2) = relatively_prime_aux2(h, fs);
            List<std::tuple<T, Integer>> new_fs = fs2;
            new_fs.insert(new_fs.begin(), std::make_tuple(f, k));
            return std::make_tuple(u, hs, new_fs);
        }
        List<T> temp1{ed::euclid_div(h, d), d};
        List<T> temp2(utils::to_unsigned_long(k), ed::euclid_div(f, d));
        List<T> temp3(utils::to_unsigned_long(k), d);
        return std::make_tuple(T(1), utils::concat(temp1, utils::concat(temp2, temp3)), fs);
    }

    template <typename T>
    std::tuple<T, List<std::tuple<T, Integer>>> relatively_prime_aux(T u, List<T> vals, List<std::tuple<T, Integer>> fs)
    {
        if (vals.size() == 0)
        {
            return std::make_tuple(u, fs);
        }
        T h = vals.at(0);
        List<T> t = utils::tail(vals);
        if (ed::is_unit(h))
        {
            return relatively_prime_aux(T(h * u), utils::tail(vals), fs);
        }
        T u2;
        List<T> hs;
        List<std::tuple<T, Integer>> fs2;
        std::tie(u2, hs, fs2) = relatively_prime_aux2(h, fs);
        return relatively_prime_aux(T(u2 * u), utils::concat(hs, t), fs2);
    }

    template <typename T>
    std::tuple<T, List<std::tuple<T, Integer>>> relatively_prime_factors(T a, T b)
    {
        return relatively_prime_aux(T(1), List<T>{a, b}, List<std::tuple<T, Integer>>{});
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

    Pair<Integer> mul(Pair<Integer> p1, Pair<Integer> p2, Integer n, Integer r, Integer s)
    {
        Integer a, b, c, d;
        std::tie(a, b) = p1;
        std::tie(c, d) = p2;
        Integer x = a * c;
        Integer y = a * d + b * c;
        Integer z = b * d;
        Integer a2 = y - x * r;
        Integer b2 = z - x * s;
        Integer a3 = utils::mod(a2, n);
        Integer b3 = utils::mod(b2, n);
        return std::make_tuple(a3, b3);
    }

    Pair<Integer> pow(Pair<Integer> x, Integer m, Integer n, Integer r, Integer s)
    {
        if (m <= 0)
        {
            return std::make_tuple(0, 1);
        }
        if (ring::odd(m))
        {
            return mul(x, pow(x, m - 1, n, r, s), n, r, s);
        }
        Pair<Integer> y = pow(x, utils::div(m, 2), n, r, s);
        return mul(y, y, n, r, s);
    };

    StepComp<Integer> root_mod(Integer n, Integer a)
    {
        if (utils::mod(a, n) == -1)
        {
            return root_of_negative_one(n);
        }

        Integer b = utils::randint(0, n - 1);
        Integer r = utils::mod(2 * b, n);
        Integer s = utils::mod(b * b - a, n);
        Integer c, d;
        std::tie(c, d) = pow(std::make_tuple(1_mpz, 0_mpz), utils::div(n - 1, 2), n, r, s);

        auto res = [=]()
        {
            std::optional<Integer> imod = ed::inv_mod(n, c);
            if (imod.has_value())
            {
                Integer c2 = imod.value();
                Integer t = (1 - d) * c2;
                Integer t1 = utils::mod(t + b, n);
                if (utils::mod(t1 * t1 - a, n) == 0)
                {
                    return StepComp<Integer>(t1);
                }
            }
            return root_mod(n, a);
        };

        return StepComp<Integer>(res);
    }

    StepComp<Maybe<ZOmega>> dioph_int_assoc_prime(Integer n)
    {
        if (n < 0)
        {
            return dioph_int_assoc_prime(-n);
        }
        if (n == 0)
        {
            return StepComp(Maybe<ZOmega>(0));
        }
        if (n == 2)
        {
            return StepComp(Maybe<ZOmega>(ring::rootTwo<ZOmega>()));
        }
        if (utils::mod(n, 4) == 1)
        {
            StepComp<Integer> s = root_of_negative_one(n);
            auto g = [=](Integer h) -> StepComp<Maybe<ZOmega>>
            {
                ZOmega t = ed::euclid_gcd<ZOmega>(
                    ring::fromInteger<ZOmega>(h) + ring::i<ZOmega>(), ring::fromInteger<ZOmega>(n));
                assert(t.adj() * t == ring::fromInteger<ZOmega>(n)); // Make sure solution is correct.
                return StepComp<Maybe<ZOmega>>(Maybe<ZOmega>(t));
            };
            return sc::bind<Integer, Maybe<ZOmega>>(s, g);
        }
        if (utils::mod(n, 8) == 3)
        {
            StepComp<Integer> s = root_mod(n, -2);
            auto g = [=](Integer h) -> StepComp<Maybe<ZOmega>>
            {
                ZOmega t = ed::euclid_gcd<ZOmega>(
                    ring::fromInteger<ZOmega>(h) + ring::i<ZOmega>() * ring::rootTwo<ZOmega>(),
                    ring::fromInteger<ZOmega>(n));
                assert(t.adj() * t == ring::fromInteger<ZOmega>(n)); // Make sure solution is correct.
                return StepComp<Maybe<ZOmega>>(Maybe<ZOmega>(t));
            };
            return sc::bind<Integer, Maybe<ZOmega>>(s, g);
        }
        if (utils::mod(n, 8) == 7)
        {
            StepComp<Integer> s = root_mod(n, 2);
            auto g = [=](Integer h __attribute__((unused))) -> StepComp<Maybe<ZOmega>>
            {
                return StepComp<Maybe<ZOmega>>(Maybe<ZOmega>());
            };
            return sc::bind<Integer, Maybe<ZOmega>>(s, g);
        }
        throw std::runtime_error("Error in dioph_int_assoc_prime");
    }

    StepComp<Maybe<ZOmega>> dioph_int_assoc_interleave(StepComp<Maybe<ZOmega>> p, StepComp<Integer> f, Integer n)
    {
        StepComp<StepComp<Maybe<ZOmega>>> p2 = p.subtask(4);
        auto g = [=](StepComp<Maybe<ZOmega>> p) -> StepComp<Maybe<ZOmega>>
        {
            if (p.is_done())
            {
                return StepComp(p.value());
            }
            StepComp<StepComp<Integer>> f2 = f.subtask(1000);
            auto g2 = [=](StepComp<Integer> f) -> StepComp<Maybe<ZOmega>>
            {
                if (f.is_done())
                {
                    Integer a = f.value();
                    int k = f.count();
                    Integer b = utils::div(n, a);
                    Integer u;
                    List<std::tuple<Integer, Integer>> facs;
                    std::tie(u, facs) = relatively_prime_factors(a, b);
                    return dioph_int_assoc_powers(facs).forward(utils::div(k, 2));
                }
                return dioph_int_assoc_interleave(p, f, n);
            };
            return sc::bind<StepComp<Integer>, Maybe<ZOmega>>(f2, g2);
        };
        return sc::bind<StepComp<Maybe<ZOmega>>, Maybe<ZOmega>>(p2, g);
    }

    StepComp<Maybe<ZOmega>> dioph_int_assoc(Integer n)
    {
        if (n < 0)
        {
            return dioph_int_assoc(-n);
        }
        if (n == 0)
        {
            return StepComp(Maybe<ZOmega>(0));
        }
        if (n == 1)
        {
            return StepComp(Maybe<ZOmega>(1));
        }
        StepComp<Maybe<ZOmega>> prime_solver = dioph_int_assoc_prime(n);
        StepComp<Integer> factor_solver = find_factor(n).speedup(30);
        return dioph_int_assoc_interleave(prime_solver, factor_solver, n);
    }

    StepComp<Maybe<ZOmega>> dioph_int_assoc_powers(List<Pair<Integer>> facs)
    {
        List<StepComp<Maybe<ZOmega>>> stepcomps;
        std::transform(facs.begin(), facs.end(), std::back_inserter(stepcomps), dioph_int_assoc_power);
        StepComp<Maybe<List<ZOmega>>> parallel = sc::parallel_list_maybe(stepcomps);
        auto g = [](Maybe<List<ZOmega>> res) -> StepComp<Maybe<ZOmega>>
        {
            if (res.has_value())
            {
                return StepComp(Maybe<ZOmega>(utils::product(res.value())));
            }
            return StepComp(Maybe<ZOmega>());
        };
        return sc::bind<Maybe<List<ZOmega>>, Maybe<ZOmega>>(parallel, g);
    }

    StepComp<Maybe<ZOmega>> dioph_int_assoc_power(Pair<Integer> p)
    {
        Integer n, k;
        std::tie(n, k) = p;
        if (ring::even(k))
        {
            ZOmega val = ring::fromInteger<ZOmega>(ring::powNonNeg<Integer>(n, utils::div(k, 2)));
            return StepComp(Maybe<ZOmega>(val));
        }
        StepComp<Maybe<ZOmega>> sc = dioph_int_assoc(n);
        auto g = [=](Maybe<ZOmega> t) -> StepComp<Maybe<ZOmega>>
        {
            if (t.has_value())
            {
                return StepComp(Maybe<ZOmega>(ring::powNonNeg(t.value(), k)));
            }
            return StepComp(Maybe<ZOmega>());
        };
        return sc::bind<Maybe<ZOmega>, Maybe<ZOmega>>(sc, g);
    }

    StepComp<Maybe<ZOmega>> dioph_zroottwo_selfassociate(ZRootTwo xi)
    {
        if (xi == 0)
        {
            return StepComp(Maybe<ZOmega>(0));
        }
        Integer a = xi.a();
        Integer b = xi.b();
        Integer n = ed::euclid_gcd(a, b);
        ZRootTwo r = ed::euclid_div(xi, ring::fromInteger<ZRootTwo>(n));
        StepComp<Maybe<ZOmega>> sc = dioph_int_assoc(n);
        auto g = [=](Maybe<ZOmega> res) -> StepComp<Maybe<ZOmega>>
        {
            if (!res.has_value())
            {
                return StepComp(Maybe<ZOmega>());
            }
            ZOmega t = res.value();
            if (ed::euclid_divides(ring::rootTwo<ZRootTwo>(), r))
            {
                return StepComp(Maybe<ZOmega>((ZOmega(1) + ring::omega<ZOmega>()) * t));
            }
            return StepComp(Maybe<ZOmega>(t));
        };
        return sc::bind<Maybe<ZOmega>, Maybe<ZOmega>>(sc, g);
    }

    StepComp<Maybe<ZOmega>> dioph_zroottwo_assoc_prime(ZRootTwo xi)
    {
        if (xi == 0)
        {
            return StepComp(Maybe<ZOmega>(0));
        }
        Integer n = ring::abs(xi.norm());
        if (utils::mod(n, 8) == 1)
        {
            StepComp<Integer> sc = root_of_negative_one(n);
            auto g = [=](Integer h) -> StepComp<Maybe<ZOmega>>
            {
                ZOmega t = ed::euclid_gcd(
                    ring::fromInteger<ZOmega>(h) + ring::i<ZOmega>(), ring::fromZRootTwo<ZOmega>(xi));
                // Make sure solution is correct.
                assert(ed::euclid_associates(t.adj() * t, ring::fromZRootTwo<ZOmega>(xi)));
                return StepComp(Maybe<ZOmega>(t));
            };
            return sc::bind<Integer, Maybe<ZOmega>>(sc, g);
        }
        if (utils::mod(n, 8) == 7)
        {
            return StepComp(Maybe<ZOmega>());
        }
        return sc::diverge<Maybe<ZOmega>>();
    }

    StepComp<Maybe<ZOmega>> dioph_zroottwo_assoc_interleave(StepComp<Maybe<ZOmega>> p, StepComp<Integer> f, ZRootTwo xi)
    {
        StepComp<StepComp<Maybe<ZOmega>>> p2 = p.subtask(4);
        auto g = [=](StepComp<Maybe<ZOmega>> p) -> StepComp<Maybe<ZOmega>>
        {
            if (p.is_done())
            {
                return StepComp(p.value());
            }
            StepComp<StepComp<Integer>> f2 = f.subtask(1000);
            auto g2 = [=](StepComp<Integer> f) -> StepComp<Maybe<ZOmega>>
            {
                if (f.is_done())
                {
                    Integer a = f.value();
                    int k = f.count();
                    ZRootTwo alpha = ed::euclid_gcd(xi, ring::fromInteger<ZRootTwo>(a));
                    ZRootTwo beta = ed::euclid_div(xi, alpha);
                    ZRootTwo u;
                    List<std::tuple<ZRootTwo, Integer>> facs;
                    std::tie(u, facs) = relatively_prime_factors(alpha, beta);
                    return dioph_zroottwo_assoc_powers(facs).forward(utils::div(k, 2));
                }
                return dioph_zroottwo_assoc_interleave(p, f, xi);
            };
            return sc::bind<StepComp<Integer>, Maybe<ZOmega>>(f2, g2);
        };
        return sc::bind<StepComp<Maybe<ZOmega>>, Maybe<ZOmega>>(p2, g);
    }

    StepComp<Maybe<ZOmega>> dioph_zroottwo_assoc(ZRootTwo xi)
    {
        if (xi == 0)
        {
            return StepComp(Maybe<ZOmega>(ZOmega(0)));
        }
        StepComp<Maybe<ZOmega>> prime_solver = dioph_zroottwo_assoc_prime(xi);
        Integer n = ring::abs(xi.norm());
        StepComp<Integer> factor_solver = find_factor(n).speedup(30);
        return dioph_zroottwo_assoc_interleave(prime_solver, factor_solver, xi);
    }

    StepComp<Maybe<ZOmega>> dioph_zroottwo_assoc_powers(List<std::tuple<ZRootTwo, Integer>> facs)
    {
        List<StepComp<Maybe<ZOmega>>> stepcomps;
        std::transform(facs.begin(), facs.end(), std::back_inserter(stepcomps), dioph_zroottwo_assoc_power);
        StepComp<Maybe<List<ZOmega>>> parallel = sc::parallel_list_maybe(stepcomps);
        auto g = [](Maybe<List<ZOmega>> res) -> StepComp<Maybe<ZOmega>>
        {
            if (!res.has_value())
            {
                return StepComp(Maybe<ZOmega>());
            }
            return StepComp(Maybe<ZOmega>(utils::product(res.value())));
        };
        return sc::bind<Maybe<List<ZOmega>>, Maybe<ZOmega>>(parallel, g);
    }

    StepComp<Maybe<ZOmega>> dioph_zroottwo_assoc_power(std::tuple<ZRootTwo, Integer> p)
    {
        ZRootTwo xi;
        Integer k;
        std::tie(xi, k) = p;
        if (ring::even(k))
        {
            return StepComp(Maybe<ZOmega>(ring::fromZRootTwo<ZOmega>(ring::powNonNeg(xi, utils::div(k, 2)))));
        }
        StepComp<Maybe<ZOmega>> sc = dioph_zroottwo_assoc(xi);
        auto g = [=](Maybe<ZOmega> t) -> StepComp<Maybe<ZOmega>>
        {
            if (!t.has_value())
            {
                return StepComp(Maybe<ZOmega>());
            }
            return StepComp(Maybe<ZOmega>(ring::powNonNeg(t.value(), k)));
        };
        return sc::bind<Maybe<ZOmega>, Maybe<ZOmega>>(sc, g);
    }
}