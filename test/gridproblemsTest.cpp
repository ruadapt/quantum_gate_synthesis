#include "../gridproblems.h"
#include <iostream>

void testLambda()
{
    std::cout << "Lambda testing:" << std::endl;

    assert(ZRootTwo(1, 1) == gridprob::lamba<ZRootTwo>());
    assert(DRootTwo(1, 1) == gridprob::lamba<DRootTwo>());
    assert(QRootTwo(1, 1) == gridprob::lamba<QRootTwo>());
    std::cout << "\tlambda tests passed" << std::endl;
}

int main()
{
    testLambda();
    std::cout << std::endl;
}