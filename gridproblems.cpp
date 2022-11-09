#include "ring.h"

namespace gridprob
{
    template <typename T>
    T lamba()
    {
        return T(1) + ring::rootTwo<T>();
    }
}