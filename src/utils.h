/** \file utils.h
 */
#pragma once
#include "types.h"
#include <gmpxx.h>
#include <algorithm>
#include <random>

/**
 * @brief Utility functions (some of which are versions of built-in Haskell functions).
 */
namespace utils
{
    /**
     * @brief Convert an Integer to an int (or throw an error if the Integer wouldn't fit).
     */
    int to_int(Integer n)
    {
        assert(n.fits_sint_p());
        return static_cast<int>(n.get_si());
    }

    /**
     * @brief Convert an Integer to an unsigned long (or throw an error if the Integer wouldn't fit).
     */
    unsigned long to_unsigned_long(Integer n)
    {
        assert(n >= 0);
        assert(n.fits_ulong_p());
        return static_cast<unsigned long>(n.get_ui());
    }

    /**
     * @brief Get the tail of a list (everything after the first element).
     * 
     * This will throw an error for an empty list.
     */
    template <typename T>
    List<T> tail(List<T> x)
    {
        if (x.size() == 0)
        {
            throw std::invalid_argument("Can't take the tail of an empty list");
        }
        return List<T>(x.begin() + 1, x.end());
    }

    /**
     * @brief Concatenate two lists into a single list.
     */
    template <typename T>
    List<T> concat(List<T> x, List<T> y)
    {
        List<T> result = x;
        result.insert(result.end(), y.begin(), y.end());
        return result;
    }

    /**
     * Flatten a list of lists into a single list.
     */
    template <typename T>
    List<T> concat(List<List<T>> list_of_lists)
    {
        List<T> result{};
        for (List<T> lst : list_of_lists)
        {
            result.insert(result.end(), lst.begin(), lst.end());
        }
        return result;
    }

    /**
     * Add an element to the front of a list.
     */
    template <typename T>
    List<T> cons(T x, List<T> lst)
    {
        return concat(List<T>{x}, lst);
    }

    template <typename A, typename B>
    List<B> map(std::function<B(A)> f, List<A> lst)
    {
        List<B> result;
        std::transform(lst.begin(), lst.end(), std::back_inserter(result), f);
        return result;
    }

    /**
     * @brief Returns the max value in a list, or the if_empty value for an empty input.
     */
    template <typename T>
    T max(List<T> lst, T if_empty)
    {
        if (lst.empty())
        {
            return if_empty;
        }
        return *std::max_element(lst.begin(), lst.end());
    }

    /**
     * This function should match the behavior of Haskell's foldr, but with an iterative
     * implementation.
     */
    template <typename A, typename B>
    B foldr(std::function<B(A, B)> f, B z, List<A> lst)
    {
        if (lst.empty())
        {
            return z;
        }
        B result = z;
        for (auto it = lst.rbegin(); it != lst.rend(); it++)
        {
            result = f(*it, result);
        }
        return result;
    }

    /**
     * This function should match the behavior of Haskell's foldl or foldl', but with an
     * iterative implementation. 
     */
    template <typename A, typename B>
    B foldl(std::function<B(B, A)> f, B z, List<A> lst)
    {
        if (lst.empty())
        {
            return z;
        }
        B result = z;
        for (auto it = lst.begin(); it != lst.end(); it++)
        {
            result = f(result, *it);
        }
        return result;
    }


    template <typename T>
    T product(List<T> list)
    {
        T result = T(1);
        for (T x : list)
        {
            result = result * x;
        }
        return result;
    }

    Integer div(Integer x, Integer y)
    {
        if (y == 0)
        {
            throw std::invalid_argument("Division by 0");
        }
        Integer q;
        mpz_fdiv_q(q.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
        return q;
    }

    int div(int x, int y)
    {
        return to_int(div(Integer(x), Integer(y)));
    }

    Integer mod(Integer x, Integer y)
    {
        Integer r;
        mpz_fdiv_r(r.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
        return r;
    }

    int mod(int x, int y)
    {
        return to_int(mod(Integer(x), Integer(y)));
    }

    /**
     * @brief Return a uniform random int in [lo, hi] (both endpoints are inclusive).
     * 
     *  If lo > hi, the two values are swapped (hi becomes lo and lo becomes hi):
     * this matches the Haskell behavior of randomR.
     */
    int randint(int lo, int hi)
    {
        if (lo > hi)
        {
            int temp = hi;
            hi = lo;
            lo = temp;
        }
        assert(lo <= hi);
        std::random_device r;
        std::default_random_engine e1(r());
        std::uniform_int_distribution<int> gen(lo, hi);
        return gen(e1);
    }

    gmp_randclass rr(gmp_randinit_default);
    bool seeded = false;

    /**
     * @brief Return a uniform random Integer in [lo, hi] (both endpoints are inclusive).
     * 
     *  If lo > hi, the two values are swapped (hi becomes lo and lo becomes hi):
     * this matches the Haskell behavior of randomR.
     */
    Integer randint(Integer lo, Integer hi)
    {
        if (lo > hi)
        {
            Integer temp = hi;
            hi = lo;
            lo = temp;
        }
        assert(lo <= hi);

        if (!seeded)
        {
            rr.seed((unsigned long)(time(NULL)));
            seeded = true;
        }

        Integer n = hi - lo + 1;
        Integer rand = rr.get_z_range(n);
        return lo + rand;
    }
}