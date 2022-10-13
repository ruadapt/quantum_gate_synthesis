#include "../ring.cpp"
#include <assert.h>
#include <iostream>
#include <memory>

template <typename T>
void testDyadic()
{
    Dyadic<T> d = Dyadic<T>(1, 5);
    Dyadic<T> d0 = Dyadic<T>(0, 9);
    Dyadic<T> dcopy = Dyadic<T>(1, 5);
    Dyadic<T> d2 = Dyadic<T>(2, 6); // 2 / 2^6 = 1 / 2^5.
    Dyadic<T> d3 = Dyadic<T>(3, 7);
    Dyadic<T> d4 = Dyadic<T>(3, 5);
    Dyadic<T> dneg = Dyadic<T>(-2, 6);

    std::cout << "Dyadic<" << typeid(T).name() << "> testing:" << std::endl;

    assert(d == dcopy); // Equal with the same constructor values.
    assert(d == d2);    // Equal with different constructor values.
    assert(!(d == d3));
    std::cout << "\tequality tests passed" << std::endl;

    Dyadic<T> d2copy = d2.copy();
    assert(d2copy == d2);
    assert(std::addressof(d2copy) != std::addressof(d2)); // Not the same object
    std::cout << "\tcopy test passed" << std::endl;

    assert("Dyadic(1, 5)" == d.toString());
    assert("Dyadic(-2, 6)" == dneg.toString());
    std::cout << "\ttoString tests passed" << std::endl;

    assert(LT == d.compare(d4));
    assert(d < d4);
    assert(d <= d4);
    assert(!(d > d4));
    assert(!(d >= d4));
    assert(EQ == d.compare(d2));
    assert(!(d > d2));
    assert(!(d < d2));
    assert(d <= d2);
    assert(d >= d2);
    assert(GT == d2.compare(d3));
    assert(d2 > d3);
    assert(d2 >= d3);
    assert(!(d2 < d3));
    assert(!(d2 <= d3));
    std::cout << "\tcomparison tests passed" << std::endl;

    assert(Dyadic<T>(7, 7) == d + d3);
    std::cout << "\tsum test passed" << std::endl;

    assert(Dyadic<T>(-1, 7) == d3 - d);
    std::cout << "\tsubtraction test passed" << std::endl;

    assert(Dyadic<T>(6, 13) == d2 * d3);
    std::cout << "\tproduct test passed" << std::endl;

    assert(Dyadic<T>(-3, 7) == -d3);
    std::cout << "\tnegation test passed" << std::endl;

    assert(1 == d.signum());
    assert(-1 == dneg.signum());
    assert(0 == d0.signum());
    std::cout << "\tsignum tests passed" << std::endl;

    assert(Dyadic<T>(1, 5) == d.abs());
    assert(Dyadic<T>(2, 6) == dneg.abs());
    std::cout << "\tabs tests passed" << std::endl;

    assert(Dyadic<T>(1, 5) == d.adj());
    assert(Dyadic<T>(-2, 6) == dneg.adj());
    std::cout << "\tadj tests passed" << std::endl;

    assert(Dyadic<T>(1, 5) == d.adj2());
    assert(Dyadic<T>(-2, 6) == dneg.adj2());
    std::cout << "\tadj2 tests passed" << std::endl;

    assert(std::make_tuple(1, 5) == d.decomposeDyadic());
    assert(std::make_tuple(-1, 5) == dneg.decomposeDyadic());
    Dyadic<T> veryReducible = Dyadic<T>(128, 9);
    assert(std::make_tuple(1, 2) == veryReducible.decomposeDyadic());
    std::cout << "\tdecomposeDyadic tests passed" << std::endl;

    assert(2 == d.integerOfDyadic(6));
    assert(-1 == dneg.integerOfDyadic(5));
    std::cout << "\tintegerOfDyadic tests passed" << std::endl;

    assert(Dyadic<T>(6, 0) == Dyadic<T>::fromInteger(6));
    std::cout << "\tfromInteger test passed" << std::endl;

    assert(Dyadic<T>(1, 5) == Dyadic<T>::fromDyadic(d));
    assert(Dyadic<T>(-2, 6) == Dyadic<T>::fromDyadic(dneg));
    std::cout << "\tfromDyadic tests passed" << std::endl;

    assert(Dyadic<T>(1, 1) == Dyadic<T>::half());
    std::cout << "\thalf test passed" << std::endl;
}

template <typename T>
void testRootTwoIntegral()
{
    // rneg2 < r < r3 < r2
    RootTwo<T> r = RootTwo<T>(1, 2);
    RootTwo<T> rcopy = RootTwo<T>(1, 2);
    RootTwo<T> r2 = RootTwo<T>(4, 9);
    RootTwo<T> r3 = RootTwo<T>(-2, 5);
    RootTwo<T> rneg1 = RootTwo<T>(-3, 1);
    RootTwo<T> rneg2 = RootTwo<T>(3, -3);

    std::cout << "RootTwo<" << typeid(T).name() << "> testing:" << std::endl;

    assert(r == rcopy);
    assert(!(r == r2));
    assert(r != r2);
    assert(!(r != rcopy));
    std::cout << "\tequality tests passed" << std::endl;

    assert("RootTwo(1, 2)" == r.toString());
    assert("RootTwo(4, 9)" == r2.toString());
    std::cout << "\ttoString tests passed" << std::endl;

    assert(rneg2 < r);
    assert(r < r3);
    assert(r3 < r2);
    assert(rneg1 < r2);
    assert(rneg2 <= r);
    assert(r <= r3);
    assert(r3 <= r2);
    assert(rneg1 <= r2);
    assert(!(rneg2 > r));
    assert(!(r > r3));
    assert(!(r3 > r2));
    assert(!(rneg1 >= r2));
    assert(!(rneg2 >= r));
    assert(!(r >= r3));
    assert(!(r3 >= r2));
    assert(!(rneg1 >= r2));
    assert(r > rneg2);
    assert(r3 > r);
    assert(r2 > r3);
    assert(r2 > rneg1);
    assert(r >= rneg2);
    assert(r3 >= r);
    assert(r2 >= r3);
    assert(r2 >= rneg1);
    assert(!(r3 < r));
    assert(!(r2 < r3));
    assert(!(r2 < rneg1));
    assert(!(r <= rneg2));
    assert(!(r3 <= r));
    assert(!(r2 <= r3));
    assert(!(r2 <= rneg1));
    std::cout << "\tcomparison tests passed" << std::endl;

    assert(RootTwo<T>(5, 2) == r + T(4));
    std::cout << "\tscalar sum test passed" << std::endl;

    assert(RootTwo<T>(-1, 2) == r - T(2));
    std::cout << "\tscalar difference test passed" << std::endl;

    assert(RootTwo<T>(12, 27) == r2 * T(3));
    std::cout << "\tscalar product test passed" << std::endl;

    assert(RootTwo<T>(5, 11) == r + r2);
    std::cout << "\tsum test passed" << std::endl;

    assert(RootTwo<T>(-3, -7) == r - r2);
    std::cout << "\tdifference test passed" << std::endl;

    assert(RootTwo<T>(40, 17) == r * r2);
    std::cout << "\tproduct test passed" << std::endl;

    assert(RootTwo<T>(-1, -2) == -r);
    std::cout << "\tnegation test passed" << std::endl;

    assert(RootTwo<T>(1, 2) == r.abs());
    assert(RootTwo<T>(-2, 5) == r3.abs());
    assert(RootTwo<T>(3, -1) == rneg1.abs());
    assert(RootTwo<T>(-3, 3) == rneg2.abs());
    std::cout << "\tabs tests passed" << std::endl;

    assert(1 == r.signum());
    assert(1 == r3.signum());
    assert(-1 == rneg1.signum());
    assert(-1 == rneg2.signum());
    std::cout << "\tsignum tests passed" << std::endl;

    assert(RootTwo<T>(T(0), 1) == RootTwo<T>::rootTwo());
    std::cout << "\trootTwo test passed" << std::endl;

    assert(RootTwo<T>(5, 0) == RootTwo<T>::fromInteger(5));
    assert(RootTwo<T>(-2, 0) == RootTwo<T>::fromInteger(-2));
    std::cout << "\tfromInteger tests passed" << std::endl;
}

template <typename T>
void testRootTwoDyadic()
{
    Dyadic<T> d0 = Dyadic<T>(0, 0);
    Dyadic<T> d1 = Dyadic<T>(1, 2);
    Dyadic<T> d2 = Dyadic<T>(2, 3);
    Dyadic<T> d3 = Dyadic<T>(-2, 5);
    Dyadic<T> d4 = Dyadic<T>(3, 5);
    Dyadic<T> d5 = Dyadic<T>(3, 7);

    RootTwo<Dyadic<T>> r0 = RootTwo<Dyadic<T>>(d0, d0);    // 0
    RootTwo<Dyadic<T>> r = RootTwo<Dyadic<T>>(d1, d2);     // 0.604
    RootTwo<Dyadic<T>> rcopy = RootTwo<Dyadic<T>>(d1, d2); // 0.604
    RootTwo<Dyadic<T>> r2 = RootTwo<Dyadic<T>>(d1, d3);    // 0.162
    RootTwo<Dyadic<T>> r3 = RootTwo<Dyadic<T>>(d5, d3);    //  -0.065
    RootTwo<Dyadic<T>> r4 = RootTwo<Dyadic<T>>(d3, d4);    // 0.070

    std::cout << "RootTwo<Dyadic<" << typeid(T).name() << ">> testing:" << std::endl;

    assert(r == rcopy);
    assert(!(r == r2));
    assert(r != r2);
    assert(!(r != rcopy));
    std::cout << "\tequality tests passed" << std::endl;

    assert("RootTwo(Dyadic(1, 2), Dyadic(2, 3))" == r.toString());
    assert("RootTwo(Dyadic(1, 2), Dyadic(-2, 5))" == r2.toString());
    std::cout << "\ttoString tests passed" << std::endl;

    // r3 < r4 < r2 < r
    assert(r3 < r4);
    assert(r4 < r2);
    assert(r2 < r);
    assert(r3 < r2);
    assert(r4 < r);
    assert(r3 <= r4);
    assert(r4 <= r2);
    assert(r2 <= r);
    assert(r3 <= r2);
    assert(r4 <= r);
    assert(!(r3 > r4));
    assert(!(r4 > r2));
    assert(!(r2 > r));
    assert(!(r3 > r2));
    assert(!(r4 > r));
    assert(!(r3 >= r4));
    assert(!(r4 >= r2));
    assert(!(r2 >= r));
    assert(!(r3 >= r2));
    assert(!(r4 >= r));
    assert(r4 > r3);
    assert(r2 > r4);
    assert(r > r2);
    assert(r2 > r3);
    assert(r > r4);
    assert(r4 >= r3);
    assert(r2 >= r4);
    assert(r >= r2);
    assert(r2 >= r3);
    assert(r >= r4);
    assert(!(r4 < r3));
    assert(!(r2 < r4));
    assert(!(r < r2));
    assert(!(r2 < r3));
    assert(!(r < r4));
    assert(!(r4 <= r3));
    assert(!(r2 <= r4));
    assert(!(r <= r2));
    assert(!(r2 <= r3));
    assert(!(r <= r4));
    std::cout << "\tcomparison tests passed" << std::endl;

    RootTwo<Dyadic<T>> scalarSum = r + Dyadic<T>(T(4), T(3));
    RootTwo<Dyadic<T>> expectedScalarSum = RootTwo<Dyadic<T>>(
        Dyadic<T>(6, 3), Dyadic<T>(2, 3));
    assert(expectedScalarSum == scalarSum);
    std::cout << "\tscalar sum test passed" << std::endl;

    RootTwo<Dyadic<T>> scalarDifference = r - Dyadic<T>(T(7), T(4));
    RootTwo<Dyadic<T>> expectedScalarDifference = RootTwo<Dyadic<T>>(
        Dyadic<T>(-3, 4), Dyadic<T>(2, 3));
    assert(expectedScalarDifference == scalarDifference);
    std::cout << "\tscalar difference test passed" << std::endl;

    RootTwo<Dyadic<T>> scalarProduct = r2 * Dyadic<T>(T(3), T(2));
    RootTwo<Dyadic<T>> expectedScalarProduct = RootTwo<Dyadic<T>>(
        Dyadic<T>(3, 4), Dyadic<T>(-6, 7));
    assert(expectedScalarProduct == scalarProduct);
    std::cout << "\tscalar product test passed" << std::endl;

    RootTwo<Dyadic<T>> sum = r + r2;
    RootTwo<Dyadic<T>> expectedSum = RootTwo<Dyadic<T>>(
        Dyadic<T>(2, 2), Dyadic<T>(6, 5));
    assert(expectedSum == sum);
    std::cout << "\tsum test passed" << std::endl;

    RootTwo<Dyadic<T>> difference = r - r2;
    RootTwo<Dyadic<T>> expectedDifference = RootTwo<Dyadic<T>>(
        Dyadic<T>(0, 2), Dyadic<T>(10, 5));
    assert(expectedDifference == difference);
    std::cout << "\tdifference test passed" << std::endl;

    RootTwo<Dyadic<T>> product = r * r2;
    RootTwo<Dyadic<T>> expectedProduct = RootTwo<Dyadic<T>>(
        Dyadic<T>(8, 8), Dyadic<T>(6, 7));
    assert(expectedProduct == product);
    std::cout << "\tproduct test passed" << std::endl;

    RootTwo<Dyadic<T>> negation = -r;
    RootTwo<Dyadic<T>> expectedNegation = RootTwo<Dyadic<T>>(
        Dyadic<T>(-1, 2), Dyadic<T>(-2, 3));
    assert(expectedNegation == negation);
    std::cout << "\tnegation test passed" << std::endl;

    RootTwo<Dyadic<T>> expectedR0Abs = RootTwo<Dyadic<T>>(
        Dyadic<T>(0, 0), Dyadic<T>(0, 0)); // No change because r0 = 0.
    assert(expectedR0Abs == r0.abs());
    RootTwo<Dyadic<T>> expectedR2Abs = RootTwo<Dyadic<T>>(
        Dyadic<T>(1, 2), Dyadic<T>(-2, 5)); // No change because r2 is positive.
    assert(expectedR2Abs == r2.abs());
    RootTwo<Dyadic<T>> expectedR3Abs = RootTwo<Dyadic<T>>(
        Dyadic<T>(-3, 7), Dyadic<T>(2, 5)); // Get -r3 because r3 is negative.
    assert(expectedR3Abs == r3.abs());
    std::cout << "\tabs tests passed" << std::endl;

    assert(0 == r0.signum());
    assert(1 == r.signum());
    assert(1 == r2.signum());
    assert(1 == r4.signum());
    assert(-1 == r3.signum());
    std::cout << "\tsignum tests passed" << std::endl;

    RootTwo<Dyadic<T>> expectedHalf = RootTwo<Dyadic<T>>(
        Dyadic<T>(1, 1), Dyadic<T>(0, 0));
    assert(expectedHalf == RootTwo<Dyadic<T>>::half());
    std::cout << "\thalf test passed" << std::endl;

    RootTwo<Dyadic<T>> expectedRootTwo = RootTwo<Dyadic<T>>(
        Dyadic<T>(0, 0), Dyadic<T>(1, 0));
    assert(expectedRootTwo == RootTwo<Dyadic<T>>::rootTwo());
    std::cout << "\trootTwo test passed" << std::endl;

    RootTwo<Dyadic<T>> expectedRootHalf = RootTwo<Dyadic<T>>(
        Dyadic<T>(0, 0), Dyadic<T>(1, 1));
    assert(expectedRootHalf == RootTwo<Dyadic<T>>::rootHalf());
    std::cout << "\trootHalf test passed" << std::endl;

    RootTwo<Dyadic<T>> expectedFrom5 = RootTwo<Dyadic<T>>(
        Dyadic<T>(5, 0), Dyadic<T>(0, 0));
    assert(expectedFrom5 == RootTwo<Dyadic<T>>::fromInteger(5));
    RootTwo<Dyadic<T>> expectedFromNeg2 = RootTwo<Dyadic<T>>(
        Dyadic<T>(-2, 0), Dyadic<T>(0, 0));
    assert(expectedFromNeg2 == RootTwo<Dyadic<T>>::fromInteger(-2));
    std::cout << "\tfromInteger tests passed" << std::endl;
}

void testZ2()
{
    Z2 x = Z2(0);
    Z2 x2 = Z2(0);
    Z2 y = Z2(1);
    Z2 y2 = Z2(1);

    std::cout << "Z2 testing:" << std::endl;

    assert(Z2(0) == Z2(0));
    assert(Z2(1) == Z2(1));
    assert(!(Z2(1) == Z2(0)));
    assert(!(Z2(0) == Z2(1)));
    assert(!(Z2(0) != Z2(0)));
    assert(!(Z2(1) != Z2(1)));
    assert(Z2(1) != Z2(0));
    assert(Z2(0) != Z2(1));
    std::cout << "\tequality tests passed" << std::endl;

    assert(Z2(0).toString() == "Z2(0)");
    assert(Z2(1).toString() == "Z2(1)");
    std::cout << "\ttoString tests passed" << std::endl;

    assert(Z2(0) == Z2(0) + Z2(0));
    assert(Z2(1) == Z2(0) + Z2(1));
    assert(Z2(1) == Z2(1) + Z2(0));
    assert(Z2(0) == Z2(1) + Z2(1));
    std::cout << "\tsum tests passed" << std::endl;

    assert(Z2(0) == Z2(0) - Z2(0));
    assert(Z2(1) == Z2(0) - Z2(1));
    assert(Z2(1) == Z2(1) - Z2(0));
    assert(Z2(0) == Z2(1) - Z2(1));
    std::cout << "\tsubtraction tests passed" << std::endl;

    assert(Z2(0) == Z2(0) * Z2(0));
    assert(Z2(0) == Z2(0) * Z2(1));
    assert(Z2(0) == Z2(1) * Z2(0));
    assert(Z2(1) == Z2(1) * Z2(1));
    std::cout << "\tproduct tests passed" << std::endl;

    assert(-Z2(0) == Z2(0));
    assert(-Z2(1) == Z2(1));
    std::cout << "\tnegation tests passed" << std::endl;

    assert(Z2(0).abs() == Z2(0));
    assert(Z2(1).abs() == Z2(1));
    std::cout << "\tabs tests passed" << std::endl;

    assert(Z2(0).signum() == 1);
    assert(Z2(1).signum() == 1);
    std::cout << "\tsignum tests passed" << std::endl;

    assert(Z2(0).adj() == Z2(0));
    assert(Z2(1).adj() == Z2(1));
    std::cout << "\tadj tests passed" << std::endl;

    assert(Z2(0).adj2() == Z2(0));
    assert(Z2(1).adj2() == Z2(1));
    std::cout << "\tadj2 tests passed" << std::endl;

    assert(Z2::fromInteger(0) == Z2(0));
    assert(Z2::fromInteger(1) == Z2(1));
    assert(Z2::fromInteger(123) == Z2(1));
    assert(Z2::fromInteger(222) == Z2(0));
    std::cout << "\tfromInteger tests passed" << std::endl;
}

int main()
{
    testRootTwoIntegral<int>();
    std::cout << std::endl;
    testRootTwoIntegral<mpz_class>();
    std::cout << std::endl;
    testRootTwoDyadic<int>();
    std::cout << std::endl;
    testRootTwoDyadic<mpz_class>();
    std::cout << std::endl;
    testDyadic<int>();
    std::cout << std::endl;
    testDyadic<mpz_class>();
    std::cout << std::endl;
    testZ2();
    return 0;
}