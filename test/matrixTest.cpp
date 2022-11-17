#include "../matrix.h"
#include "../ring.h"
#include <iostream>

int main()
{
    Vector<QRootTwo, 7> v;
    v(1) = QRootTwo(5_mpq / 6, 7_mpq / 8);
    std::cout << "v(1) = " << v(1).toString() << std::endl;
    std::cout << "v(2) = " << v(2).toString() << std::endl;

    Matrix<ZRootTwo, 4, 4> m;
    m(1, 1) = ZRootTwo(4, 6);
    std::cout << "m(1, 1) = " << m(1, 1).toString() << std::endl;
    std::cout << "m(2, 2) = " << m(2, 2).toString() << std::endl;

    return 0;
}