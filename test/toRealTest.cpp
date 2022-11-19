#include "../toReal.h"

void toRealTest()
{
    std::cout << "toReal testing:" << std::endl;
    assert(1.25 == toReal<double>(5_mpq / 4));
    assert(12.0 == toReal<double>(12_mpz));
    assert(12.0 == toReal<double>(12));
    assert(1.56 == toReal<double>(1.56));
    std::cout << "\ttoReal tests passed" << std::endl;
}

int main()
{
    toRealTest();
    std::cout << std::endl;
}