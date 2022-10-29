#include "../ring.h"
#include <assert.h>
#include <iostream>

void templateError(std::string T, std::string baseType)
{
    std::cout << "Unexpected template type " << T << " for " << baseType << std::endl;
    exit(1);
}

void testUtilityFunctions()
{
    std::cout << "Utility function testing:" << std::endl;

    assert(2 == ring::intsqrt(4));
    assert(11096085937082_mpz == ring::intsqrt(123123123123123123123123123_mpz));
    assert(456456456456456456_mpz == ring::intsqrt(208352496640784928656512368224079936_mpz));
    std::cout << "\tintsqrt tests passed" << std::endl;
}

void testTypeConversions()
{
    std::cout << "Type conversion testing" << std::endl;

    double d = double(2);
    assert(2.0 == d);
    Integer i = Integer(2);
    assert(2_mpz == i);
    Rational r = Rational(2);
    assert(2_mpq == r);
    Dyadic<int> dInt = Dyadic<int>(2);
    assert(Dyadic<int>(2, 0) == dInt);
    Dyadic<Integer> dInteger = Dyadic<Integer>(2);
    assert(Dyadic<Integer>(2_mpz, 0_mpz) == dInteger);
    RootTwo<int> rInt = RootTwo<int>(2);
    assert(RootTwo<int>(2, 0) == rInt);
    RootTwo<Integer> rInteger = RootTwo<Integer>(2);
    assert(RootTwo<Integer>(2_mpz, 0_mpz) == rInteger);
    RootTwo<double> rDouble = RootTwo<double>(2);
    assert(RootTwo<double>(2.0, 0.0) == rDouble);
    RootTwo<Rational> rRational = RootTwo<Rational>(2);
    assert(RootTwo<Rational>(2_mpq, 0_mpq) == rRational);
    RootTwo<Dyadic<int>> rDyadicInt = RootTwo<Dyadic<int>>(2);
    assert(RootTwo<Dyadic<int>>(Dyadic<int>(2, 0), Dyadic<int>(0, 0)) == rDyadicInt);
    RootTwo<Dyadic<Integer>> rDyadicInteger = RootTwo<Dyadic<Integer>>(2);
    assert(RootTwo<Dyadic<Integer>>(Dyadic<Integer>(2_mpz, 0_mpz), Dyadic<Integer>(0_mpz, 0_mpz)) == rDyadicInteger);
    Complex<QRootTwo> cQRootTwo = Complex<QRootTwo>(2);
    assert(Complex<QRootTwo>(QRootTwo(2_mpq, 0_mpq), QRootTwo(0_mpq, 0_mpq)) == cQRootTwo);
    Complex<DRootTwo> cDRootTwo = Complex<DRootTwo>(2);
    assert(Complex<DRootTwo>(DRootTwo(ZDyadic(2, 0), ZDyadic(0, 0)), DRootTwo(ZDyadic(0, 0), ZDyadic(0, 0))) == cDRootTwo);
    Omega<int> oInt = Omega<int>(2);
    assert(Omega<int>(0, 0, 0, 2) == oInt);
    Omega<Integer> oInteger = Omega<Integer>(2);
    assert(Omega<Integer>(0_mpz, 0_mpz, 0_mpz, 2_mpz) == oInteger);
    Omega<Rational> oRational = Omega<Rational>(2);
    assert(Omega<Rational>(0_mpq, 0_mpq, 0_mpq, 2_mpq) == oRational);
    Omega<ZDyadic> oDyadicInteger = Omega<ZDyadic>(2);
    assert(Omega<ZDyadic>(ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(2, 0)) == oDyadicInteger);
    std::cout << "\tconversion from int tests passed" << std::endl;
}

void testHalfRing()
{
    std::cout << "HalfRing testing:" << std::endl;

    assert(0.5 == ring::half<double>());
    assert(1_mpq / 2 == ring::half<Rational>());
    assert(Dyadic<int>(1, 1) == ring::half<Dyadic<int>>());
    assert(Dyadic<Integer>(1, 1) == ring::half<Dyadic<Integer>>());
    assert(RootTwo<Rational>(1_mpq / 2, 0_mpq) == ring::half<RootTwo<Rational>>());
    assert(Complex<double>(0.5, 0) == ring::half<Complex<double>>());
    assert(Omega<double>(0, 0, 0, 0.5) == ring::half<Omega<double>>());
    assert(Omega<Rational>(0_mpq, 0_mpq, 0_mpq, 1_mpq / 2) == ring::half<Omega<Rational>>());
    assert(Omega<ZDyadic>(ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(1, 1)) == ring::half<Omega<ZDyadic>>());
    std::cout << "\thalf tests passed" << std::endl;

    assert(0.625 == ring::fromDyadic<double>(Dyadic<int>(5, 3)));
    assert(40 == ring::fromDyadic<double>(Dyadic<int>(5, -3)));
    assert(5_mpq / 8 == ring::fromDyadic<Rational>(Dyadic<int>(5, 3)));
    assert(40_mpq == ring::fromDyadic<Rational>(Dyadic<int>(5, -3)));
    assert(RootTwo<Rational>(5_mpq / 8, 0_mpq) == ring::fromDyadic<RootTwo<Rational>>(Dyadic<int>(5, 3)));
    assert(RootTwo<Rational>(40_mpq, 0_mpq) == ring::fromDyadic<RootTwo<Rational>>(Dyadic<int>(5, -3)));
    assert(Complex<double>(0.625, 0) == ring::fromDyadic<Complex<double>>(Dyadic<int>(5, 3)));
    assert(Complex<double>(40, 0) == ring::fromDyadic<Complex<double>>(Dyadic<int>(5, -3)));
    assert(Omega<double>(0, 0, 0, 0.625) == ring::fromDyadic<Omega<double>>(Dyadic<int>(5, 3)));
    assert(Omega<double>(0, 0, 0, 40) == ring::fromDyadic<Omega<double>>(Dyadic<int>(5, -3)));
    assert(Omega<Rational>(0_mpq, 0_mpq, 0_mpq, 5_mpq / 8) ==
           ring::fromDyadic<Omega<Rational>>(Dyadic<int>(5, 3)));
    assert(Omega<Rational>(0_mpq, 0_mpq, 0_mpq, 40_mpq) ==
           ring::fromDyadic<Omega<Rational>>(Dyadic<int>(5, -3)));
    assert(Omega<ZDyadic>(ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(5, 3)) ==
           ring::fromDyadic<Omega<ZDyadic>>(Dyadic<int>(5, 3)));
    assert(Omega<ZDyadic>(ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(5, -3)) ==
           ring::fromDyadic<Omega<ZDyadic>>(Dyadic<int>(5, -3)));
    std::cout << "\tfromDyadic for Dyadic<int> passed" << std::endl;

    assert(0.625 == ring::fromDyadic<double>(ZDyadic(5, 3)));
    assert(40 == ring::fromDyadic<double>(ZDyadic(5, -3)));
    assert(5_mpq / 8 == ring::fromDyadic<Rational>(ZDyadic(5, 3)));
    assert(40_mpq == ring::fromDyadic<Rational>(ZDyadic(5, -3)));
    assert(RootTwo<Rational>(5_mpq / 8, 0_mpq) == ring::fromDyadic<RootTwo<Rational>>(ZDyadic(5, 3)));
    assert(RootTwo<Rational>(40_mpq, 0_mpq) == ring::fromDyadic<RootTwo<Rational>>(ZDyadic(5, -3)));
    assert(Complex<double>(0.625, 0) == ring::fromDyadic<Complex<double>>(ZDyadic(5, 3)));
    assert(Complex<double>(40, 0) == ring::fromDyadic<Complex<double>>(ZDyadic(5, -3)));
    assert(Omega<double>(0, 0, 0, 0.625) == ring::fromDyadic<Omega<double>>(ZDyadic(5, 3)));
    assert(Omega<double>(0, 0, 0, 40) == ring::fromDyadic<Omega<double>>(ZDyadic(5, -3)));
    assert(Omega<Rational>(0_mpq, 0_mpq, 0_mpq, 5_mpq / 8) ==
           ring::fromDyadic<Omega<Rational>>(ZDyadic(5, 3)));
    assert(Omega<Rational>(0_mpq, 0_mpq, 0_mpq, 40_mpq) ==
           ring::fromDyadic<Omega<Rational>>(ZDyadic(5, -3)));
    assert(Omega<ZDyadic>(ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(5, 3)) ==
           ring::fromDyadic<Omega<ZDyadic>>(ZDyadic(5, 3)));
    assert(Omega<ZDyadic>(ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(5, -3)) ==
           ring::fromDyadic<Omega<ZDyadic>>(ZDyadic(5, -3)));
    std::cout << "\tfromDyadic for Dyadic<Integer> (ZDyadic) passed" << std::endl;
}

void testRootTwoRing()
{
    std::cout << "RootTwoRing testing:" << std::endl;

    assert(1.4142135623730951 == ring::rootTwo<double>());
    assert(Complex<double>(1.4142135623730951, 0) == ring::rootTwo<Complex<double>>());
    assert(RootTwo<Integer>(0, 1) == ring::rootTwo<RootTwo<Integer>>());
    assert(RootTwo<Rational>(0_mpq, 1_mpq) == ring::rootTwo<RootTwo<Rational>>());
    assert(Omega<Integer>(-1, 0, 1, 0) == ring::rootTwo<Omega<Integer>>());
    assert(Omega<ZDyadic>(-1, 0, 1, 0) == ring::rootTwo<Omega<ZDyadic>>());
    assert(Omega<Rational>(-1, 0, 1, 0) == ring::rootTwo<Omega<Rational>>());
    std::cout << "\trootTwo tests passed" << std::endl;

    ZRootTwo z = ZRootTwo(10, -7);
    assert(0.10050506338833465 == ring::fromZRootTwo<double>(z));
    assert(Complex<double>(0.10050506338833465, 0) == ring::fromZRootTwo<Complex<double>>(z));
    assert(RootTwo<Integer>(10, -7) == ring::fromZRootTwo<RootTwo<Integer>>(z));
    assert(RootTwo<Rational>(10_mpq, -7_mpq) == ring::fromZRootTwo<RootTwo<Rational>>(z));
    assert(Omega<Integer>(7, 0, -7, 10) == ring::fromZRootTwo<Omega<Integer>>(z));
    assert(Omega<ZDyadic>(7, 0, -7, 10) == ring::fromZRootTwo<Omega<ZDyadic>>(z));
    assert(Omega<Rational>(7, 0, -7, 10) == ring::fromZRootTwo<Omega<Rational>>(z));
    std::cout << "\tfromZRootTwo tests passed" << std::endl;
}

void testRootHalfRing()
{
    std::cout << "RootHalfRing testing:" << std::endl;

    assert(0.7071067811865476 == ring::rootHalf<double>());
    assert(Complex<double>(0.7071067811865476, 0) == ring::rootHalf<Complex<double>>());
    assert(RootTwo<Rational>(0, 1_mpq / 2) == ring::rootHalf<RootTwo<Rational>>());
    assert(Omega<ZDyadic>(ZDyadic(-1, 1), 0, ZDyadic(1, 1), 0) == ring::rootHalf<Omega<ZDyadic>>());
    assert(Omega<Rational>(-1_mpq / 2, 0, 1_mpq / 2, 0) == ring::rootHalf<Omega<Rational>>());
    std::cout << "\trootHalf tests passed" << std::endl;

    DRootTwo d = DRootTwo(Dyadic<Integer>(5, 2), Dyadic<Integer>(-3, 5));
    assert(1.1174174785275224 == ring::fromDRootTwo<double>(d));
    assert(Complex<double>(1.1174174785275224, 0) == ring::fromDRootTwo<Complex<double>>(d));
    assert(RootTwo<Rational>(5_mpq / 4, -3_mpq / 32) == ring::fromDRootTwo<RootTwo<Rational>>(d));
    assert(Omega<ZDyadic>(ZDyadic(3, 5), 0, ZDyadic(-3, 5), ZDyadic(5, 2)) == ring::fromDRootTwo<Omega<ZDyadic>>(d));
    assert(Omega<Rational>(3_mpq / 32, 0, -3_mpq / 32, 5_mpq / 4) == ring::fromDRootTwo<Omega<Rational>>(d));
    std::cout << "\tfromDRootTwo tests passed" << std::endl;
}

void testComplexRing()
{
    std::cout << "ComplexRing testing:" << std::endl;

    assert(Complex<int>(0, 1) == ring::i<Complex<int>>());
    assert(Complex<Integer>(0, 1) == ring::i<Complex<Integer>>());
    assert(Complex<Complex<int>>(Complex<int>(0, 0), Complex<int>(1, 0)) == ring::i<Complex<Complex<int>>>());
    assert(RootTwo<Complex<int>>(Complex<int>(0, 1), Complex<int>(0, 0)) == ring::i<RootTwo<Complex<int>>>());
    assert(RootTwo<Complex<Integer>>(Complex<Integer>(0, 1), Complex<Integer>(0, 0)) == ring::i<RootTwo<Complex<Integer>>>());
    assert(Omega<Integer>(0, 1, 0, 0) == ring::i<Omega<Integer>>());
    assert(Omega<ZDyadic>(0, 1, 0, 0) == ring::i<Omega<ZDyadic>>());
    assert(Omega<Rational>(0, 1, 0, 0) == ring::i<Omega<Rational>>());
    std::cout << "\ti tests passed" << std::endl;
}

void testOmegaRing()
{
    std::cout << "OmegaRing testing:" << std::endl;

    assert(CDouble(0.7071067811865476, 0.7071067811865476) == ring::omega<CDouble>());
    assert(DRComplex(DRootTwo(0, ZDyadic(1, 1)), DRootTwo(0, ZDyadic(1, 1))) ==
           ring::omega<Complex<DRootTwo>>());
    assert(QRComplex(QRootTwo(0, 1_mpq / 2), QRootTwo(0, 1_mpq / 2)) ==
           ring::omega<Complex<QRootTwo>>());

    assert(RootTwo<DComplex>(DComplex(0, 0), DComplex(ZDyadic(1, 1), ZDyadic(1, 1))) ==
           ring::omega<RootTwo<DComplex>>());

    assert(ZOmega(0, 0, 1, 0) == ring::omega<ZOmega>());
    assert(DOmega(0, 0, 1, 0) == ring::omega<DOmega>());
    assert(QOmega(0, 0, 1, 0) == ring::omega<QOmega>());

    std::cout << "\tomega tests passed" << std::endl;
}

void testNormedRing()
{
    std::cout << "NormedRing testing:" << std::endl;

    assert(10_mpz == ring::norm<int>(10));
    assert(10_mpz == ring::norm<Integer>(10));
    assert(85_mpz == ring::norm<Complex<int>>(Complex<int>(6, -7)));
    assert(85_mpz == ring::norm<Complex<Integer>>(Complex<Integer>(6, -7)));
    assert(-62_mpz == ring::norm<RootTwo<int>>(RootTwo<int>(6, -7)));
    assert(-62_mpz == ring::norm<RootTwo<Integer>>(RootTwo<Integer>(6, -7)));
    // norm = (3^2 + 4^2)^2 - 2 * (5^2 + 2^2)^2 = 35^2 - 2 * 29^2 = -1057.
    assert(-1057_mpz == ring::norm<RootTwo<Complex<int>>>(
                            RootTwo<Complex<int>>(Complex<int>(3, 4), Complex<int>(5, 2))));
    // norm = (3^2 - 2 * 4^2)^2 + (5^2 - 2 * 2^2)^2 = (-23)^2 + (17)^2 = 818.
    assert(818_mpq == ring::norm<Complex<RootTwo<int>>>(
                          Complex<RootTwo<int>>(RootTwo<int>(3, 4), RootTwo<int>(5, 2))));
    assert(6562_mpq == ring::norm<Omega<int>>(Omega<int>(7, 3, -2, 6)));
    assert(6562_mpq == ring::norm<Omega<Integer>>(Omega<Integer>(7, 3, -2, 6)));
    std::cout << "\tnorm tests passed" << std::endl;
}

void testAdjoint()
{
    std::cout << "Adjoint testing:" << std::endl;

    assert(7 == ring::adj<int>(7));
    assert(123123123123123123_mpz == ring::adj<Integer>(123123123123123123_mpz));
    assert(0.123 == ring::adj<double>(0.123));
    assert(123_mpq / 456 == ring::adj<Rational>(123_mpq / 456));
    assert(Z2(0) == ring::adj<Z2>(0));
    assert(Z2(1) == ring::adj<Z2>(1));
    assert(Dyadic<int>(2, 7) == ring::adj<Dyadic<int>>(Dyadic<int>(2, 7)));
    assert(Dyadic<Integer>(2, 7) == ring::adj<Dyadic<Integer>>(Dyadic<Integer>(2, 7)));
    assert(Complex<int>(3, -4) == ring::adj<Complex<int>>(Complex<int>(3, 4)));
    assert(Complex<Integer>(3, -4) == ring::adj<Complex<Integer>>(Complex<Integer>(3, 4)));
    assert(RootTwo<int>(4, 7) == ring::adj<RootTwo<int>>(RootTwo<int>(4, 7)));
    assert(RootTwo<Integer>(4, 7) == ring::adj<RootTwo<Integer>>(RootTwo<Integer>(4, 7)));
    assert(
        RootTwo<Complex<int>>(Complex<int>(4, -5), Complex<int>(6, -7)) ==
        ring::adj<RootTwo<Complex<int>>>(RootTwo<Complex<int>>(Complex<int>(4, 5), Complex<int>(6, 7))));
    assert(
        RootTwo<Complex<Integer>>(Complex<Integer>(4, -5), Complex<Integer>(6, -7)) ==
        ring::adj<RootTwo<Complex<Integer>>>(RootTwo<Complex<Integer>>(Complex<Integer>(4, 5), Complex<Integer>(6, 7))));
    assert(
        Complex<RootTwo<Dyadic<int>>>(
            RootTwo<Dyadic<int>>(Dyadic<int>(5, 7), Dyadic<int>(9, 4)),
            RootTwo<Dyadic<int>>(Dyadic<int>(-7, 9), Dyadic<int>(-11, 4))) ==
        ring::adj<Complex<RootTwo<Dyadic<int>>>>(Complex<RootTwo<Dyadic<int>>>(
            RootTwo<Dyadic<int>>(Dyadic<int>(5, 7), Dyadic<int>(9, 4)),
            RootTwo<Dyadic<int>>(Dyadic<int>(7, 9), Dyadic<int>(11, 4)))));
    assert(ZOmega(-3, -2, -1, 4) == ring::adj<ZOmega>(ZOmega(1, 2, 3, 4)));
    assert(DOmega(ZDyadic(-5, 6), ZDyadic(-3, 4), ZDyadic(-1, 2), ZDyadic(7, 8)) ==
           ring::adj<DOmega>(DOmega(ZDyadic(1, 2), ZDyadic(3, 4), ZDyadic(5, 6), ZDyadic(7, 8))));
    assert(QOmega(-5_mpq / 6, -3_mpq / 4, -1_mpq / 2, 7_mpq / 8) ==
           ring::adj<QOmega>(QOmega(1_mpq / 2, 3_mpq / 4, 5_mpq / 6, 7_mpq / 8)));
    std::cout << "\tadj tests passed" << std::endl;
}

void testAdjoint2()
{
    std::cout << "Adjoint2 testing:" << std::endl;

    assert(7 == ring::adj2<int>(7));
    assert(123123123123123123_mpz == ring::adj2<Integer>(123123123123123123_mpz));
    assert(123_mpq / 456 == ring::adj2<Rational>(123_mpq / 456));
    assert(Z2(0) == ring::adj2<Z2>(0));
    assert(Z2(1) == ring::adj2<Z2>(1));
    assert(Dyadic<int>(2, 7) == ring::adj2<Dyadic<int>>(Dyadic<int>(2, 7)));
    assert(Dyadic<Integer>(2, 7) == ring::adj2<Dyadic<Integer>>(Dyadic<Integer>(2, 7)));
    assert(Complex<int>(3, 4) == ring::adj2<Complex<int>>(Complex<int>(3, 4)));
    assert(Complex<Integer>(3, 4) == ring::adj2<Complex<Integer>>(Complex<Integer>(3, 4)));
    assert(RootTwo<int>(4, -7) == ring::adj2<RootTwo<int>>(RootTwo<int>(4, 7)));
    assert(RootTwo<Integer>(4, -7) == ring::adj2<RootTwo<Integer>>(RootTwo<Integer>(4, 7)));
    assert(
        RootTwo<Complex<int>>(Complex<int>(4, 5), Complex<int>(-6, -7)) ==
        ring::adj2<RootTwo<Complex<int>>>(RootTwo<Complex<int>>(Complex<int>(4, 5), Complex<int>(6, 7))));
    assert(
        RootTwo<Complex<Integer>>(Complex<Integer>(4, 5), Complex<Integer>(-6, -7)) ==
        ring::adj2<RootTwo<Complex<Integer>>>(RootTwo<Complex<Integer>>(Complex<Integer>(4, 5), Complex<Integer>(6, 7))));
    assert(
        Complex<RootTwo<Dyadic<int>>>(
            RootTwo<Dyadic<int>>(Dyadic<int>(5, 7), Dyadic<int>(-9, 4)),
            RootTwo<Dyadic<int>>(Dyadic<int>(7, 9), Dyadic<int>(-11, 4))) ==
        ring::adj2<Complex<RootTwo<Dyadic<int>>>>(Complex<RootTwo<Dyadic<int>>>(
            RootTwo<Dyadic<int>>(Dyadic<int>(5, 7), Dyadic<int>(9, 4)),
            RootTwo<Dyadic<int>>(Dyadic<int>(7, 9), Dyadic<int>(11, 4)))));
    assert(Omega<int>(-1, 2, -3, 4) == ring::adj2<Omega<int>>(Omega<int>(1, 2, 3, 4)));
    assert(ZOmega(-1, 2, -3, 4) == ring::adj2<ZOmega>(ZOmega(1, 2, 3, 4)));
    assert(DOmega(ZDyadic(-1, 2), ZDyadic(3, 4), ZDyadic(-5, 6), ZDyadic(7, 8)) ==
           ring::adj2<DOmega>(DOmega(ZDyadic(1, 2), ZDyadic(3, 4), ZDyadic(5, 6), ZDyadic(7, 8))));
    assert(QOmega(-1_mpq / 2, 3_mpq / 4, -5_mpq / 6, 7_mpq / 8) ==
           ring::adj2<QOmega>(QOmega(1_mpq / 2, 3_mpq / 4, 5_mpq / 6, 7_mpq / 8)));
    std::cout << "\tadj2 tests passed" << std::endl;
}

void testFloor()
{
    std::cout << "Floor testing:" << std::endl;

    // double
    assert(5 == ring::floor_of(5.2));
    assert(4 == ring::floor_of(4));
    assert(-3 == ring::floor_of(-2.5));
    // Integer
    assert(123123123123123123123123_mpq == ring::floor_of(Integer(123123123123123123123123_mpz)));
    assert(-22 == ring::floor_of(Integer(-22_mpq)));
    // Rational
    assert(4 == ring::floor_of(Rational(123_mpq / 25)));
    assert(12 == ring::floor_of(Rational(144_mpq / 12)));
    assert(-13 == ring::floor_of(Rational(-145_mpq / 12)));
    // QRootTwo
    assert(3 == ring::floor_of(QRootTwo(123_mpq / 22, -45_mpq / 27)));
    assert(12 == ring::floor_of(QRootTwo(144_mpq / 12, 0)));
    assert(-1398 == ring::floor_of(QRootTwo(12345671234567_mpq / 786876876876_mpq, -999)));
    std::cout << "\tfloor_of tests passed" << std::endl;

    // double
    assert(6 == ring::ceiling_of(5.2));
    assert(4 == ring::ceiling_of(4));
    assert(-2 == ring::ceiling_of(-2.5));
    // Integer
    assert(123123123123123123123123_mpz == ring::ceiling_of(Integer(123123123123123123123123_mpz)));
    assert(-22 == ring::ceiling_of(Integer(-22_mpz)));
    // Rational
    assert(5 == ring::ceiling_of(Rational(123_mpq / 25)));
    assert(12 == ring::ceiling_of(Rational(144_mpq / 12)));
    assert(-12 == ring::ceiling_of(Rational(-145_mpq / 12)));
    // QRootTwo
    assert(4 == ring::ceiling_of(QRootTwo(123_mpq / 22, -45_mpq / 27)));
    assert(12 == ring::ceiling_of(QRootTwo(144_mpq / 12, 0)));
    assert(-1397 == ring::ceiling_of(QRootTwo(12345671234567_mpq / 786876876876_mpz, -999)));
    std::cout << "\tceiling_of tests passed" << std::endl;
}

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

    Dyadic<T> dnegDenom = Dyadic<T>(3, -5);
    assert(std::make_tuple(96, 0) == dnegDenom.decomposeDyadic());
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
    std::cout << "\tfromDyadic<Dyadic<T>> tests passed" << std::endl;

    assert(0.25 == ring::fromDyadic<double>(Dyadic<T>(1, 2)));
    std::cout << "\tfromDyadic<double> test passed" << std::endl;

    assert(3_mpq / 128 == ring::fromDyadic<Rational>(Dyadic<T>(3, 7)));
    std::cout << "\tfromDyadic<Rational> test passed" << std::endl;

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

    assert(RootTwo<T>(5, 2) == r + 4);
    std::cout << "\tscalar sum test passed" << std::endl;

    assert(RootTwo<T>(-1, 2) == r - 2);
    std::cout << "\tscalar difference test passed" << std::endl;

    assert(RootTwo<T>(12, 27) == r2 * 3);
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

    assert(RootTwo<T>(0, 1) == RootTwo<T>::rootTwo());
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

void testRootTwoRational()
{
    QRootTwo r0 = QRootTwo(0_mpq, 0_mpq);              // 0
    QRootTwo r1 = QRootTwo(2_mpq / 4, 9_mpq / 18);     // 1.207
    QRootTwo r1equal = QRootTwo(1_mpq / 2, 1_mpq / 2); // 1.207
    QRootTwo r2 = QRootTwo(-7_mpq / 13, 33_mpq / 20);  // 1.795
    QRootTwo r3 = QRootTwo(9_mpq / 17, -100_mpq / 3);  // -46.611
    QRootTwo r4 = QRootTwo(5_mpq, -7_mpq);             // -4.899

    std::cout << "QRootTwo testing:" << std::endl;

    // Make sure that we're converting the rational numbers to the standard,
    // canonical form.
    QRootTwo rNonCanonical = QRootTwo(2_mpq / 4, -3_mpq / 9);
    assert(rNonCanonical.a.get_num() == 1);
    assert(rNonCanonical.a.get_den() == 2);
    assert(rNonCanonical.b.get_num() == -1);
    assert(rNonCanonical.b.get_den() == 3);
    std::cout << "\tcanonicalization test passed" << std::endl;

    assert(r1 == r1equal);
    assert(r1 != r0);
    assert(r1 != r2);
    assert(r1 != r3);
    assert(r2 == r2);
    assert(r2 != r3);
    assert(!(r1 != r1equal));
    assert(!(r1 == r0));
    assert(!(r1 == r2));
    assert(!(r1 == r3));
    assert(!(r2 != r2));
    assert(!(r2 == r3));
    std::cout << "\tequality tests passed" << std::endl;

    assert("RootTwo(0, 0)" == r0.toString());
    assert("RootTwo(-7/13, 33/20)" == r2.toString());
    assert("RootTwo(9/17, -100/3)" == r3.toString());
    assert("RootTwo(5, -7)" == r4.toString());
    std::cout << "\ttoString tests passed" << std::endl;

    // // r3 < r4 < r0 < r1 < r2
    assert(r3 < r4);
    assert(r4 < r0);
    assert(r0 < r1);
    assert(r1 < r2);
    assert(r4 < r1);
    assert(r0 < r2);
    assert(r3 <= r4);
    assert(r4 <= r0);
    assert(r0 <= r1);
    assert(r1 <= r2);
    assert(r4 <= r1);
    assert(r0 <= r2);
    assert(!(r3 > r4));
    assert(!(r4 > r0));
    assert(!(r0 > r1));
    assert(!(r1 > r2));
    assert(!(r4 > r1));
    assert(!(r0 > r2));
    assert(!(r3 >= r4));
    assert(!(r4 >= r0));
    assert(!(r0 >= r1));
    assert(!(r1 >= r2));
    assert(!(r4 >= r1));
    assert(!(r0 >= r2));
    assert(r4 > r3);
    assert(r0 > r4);
    assert(r1 > r0);
    assert(r2 > r1);
    assert(r1 > r4);
    assert(r2 > r0);
    assert(r4 >= r3);
    assert(r0 >= r4);
    assert(r1 >= r0);
    assert(r2 >= r1);
    assert(r1 >= r4);
    assert(r2 >= r0);
    assert(!(r4 < r3));
    assert(!(r0 < r4));
    assert(!(r1 < r0));
    assert(!(r2 < r1));
    assert(!(r1 < r4));
    assert(!(r2 < r0));
    assert(!(r4 <= r3));
    assert(!(r0 <= r4));
    assert(!(r1 <= r0));
    assert(!(r2 <= r1));
    assert(!(r1 <= r4));
    assert(!(r2 <= r0));
    std::cout << "\tcomparison tests passed" << std::endl;

    QRootTwo sum = r1 + r2;
    QRootTwo expectedSum = QRootTwo(-1_mpq / 26, 43_mpq / 20);
    assert(expectedSum == sum);
    std::cout << "\tsum test passed" << std::endl;

    QRootTwo difference = r2 - r3;
    QRootTwo expectedDifference = QRootTwo(-236_mpq / 221, 2099_mpq / 60);
    assert(expectedDifference == difference);
    std::cout << "\tdifference test passed" << std::endl;

    QRootTwo product = r2 * r3;
    QRootTwo expectedProduct = QRootTwo(-24373_mpq / 221, 249583_mpq / 13260);
    assert(expectedProduct == product);
    std::cout << "\tproduct test passed" << std::endl;

    QRootTwo negation = -r3;
    QRootTwo expectedNegation = QRootTwo(-9_mpq / 17, 100_mpq / 3);
    assert(expectedNegation == negation);
    std::cout << "\tnegation test passed" << std::endl;

    QRootTwo expectedR0Abs = QRootTwo(0_mpq, 0_mpq);
    assert(QRootTwo(0_mpq, 0_mpq) == r0.abs());
    assert(QRootTwo(1_mpq / 2, 1_mpq / 2) == r1.abs());
    assert(QRootTwo(-9_mpq / 17, 100_mpq / 3) == r3.abs());
    std::cout << "\tabs tests passed" << std::endl;

    assert(0 == r0.signum());
    assert(1 == r1.signum());
    assert(1 == r2.signum());
    assert(-1 == r3.signum());
    assert(-1 == r4.signum());
    std::cout << "\tsignum tests passed" << std::endl;

    assert(QRootTwo(1_mpq / 2, 0_mpq) == QRootTwo::half());
    std::cout << "\thalf test passed" << std::endl;

    assert(QRootTwo(0_mpq, 1_mpq) == QRootTwo::rootTwo());
    std::cout << "\trootTwo test passed" << std::endl;

    assert(QRootTwo(0_mpq, 1_mpq / 2) == QRootTwo::rootHalf());
    std::cout << "\trootHalf test passed" << std::endl;

    assert(QRootTwo(5_mpq, 0_mpq) == QRootTwo::fromInteger(5));
    assert(QRootTwo(-2_mpq, 0_mpq) == QRootTwo::fromInteger(-2));
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

template <typename T>
void testComplex()
{
    Complex<T> c1 = Complex<T>(1, 2);
    Complex<T> c1copy = Complex<T>(1, 2);
    Complex<T> c2 = Complex<T>(4, 9);
    Complex<T> c3 = Complex<T>(-2, 5);
    Complex<T> c4 = Complex<T>(-3, -1);

    std::cout << "Complex<" << typeid(T).name() << "> testing:" << std::endl;

    assert(c1 == c1copy);
    assert(!(c1 == c2));
    assert(c1 != c2);
    assert(!(c1 != c1copy));
    std::cout << "\tequality tests passed" << std::endl;

    if constexpr (std::is_same<T, int>::value || std::is_same<T, Integer>::value || std::is_same<T, Rational>::value)
    {
        assert("Complex(1, 2)" == c1.toString());
        assert("Complex(4, 9)" == c2.toString());
        std::cout << "\ttoString tests passed" << std::endl;
    }
    else if constexpr (std::is_same<T, double>::value)
    {
        assert("Complex(1.000000, 2.000000)" == c1.toString());
        assert("Complex(4.000000, 9.000000)" == c2.toString());
        std::cout << "\ttoString tests passed" << std::endl;
    }
    else if constexpr (std::is_same<T, ZDyadic>::value)
    {
        assert("Complex(Dyadic(1, 0), Dyadic(2, 0))" == c1.toString());
        assert("Complex(Dyadic(4, 0), Dyadic(9, 0))" == c2.toString());
        std::cout << "\ttoString tests passed" << std::endl;
    }
    else if constexpr (std::is_same<T, DRootTwo>::value)
    {
        assert("Complex(RootTwo(Dyadic(1, 0), Dyadic(0, 0)), RootTwo(Dyadic(2, 0), Dyadic(0, 0)))" == c1.toString());
        assert("Complex(RootTwo(Dyadic(4, 0), Dyadic(0, 0)), RootTwo(Dyadic(9, 0), Dyadic(0, 0)))" == c2.toString());
        std::cout << "\ttoString tests passed" << std::endl;
    }
    else if constexpr (std::is_same<T, QRootTwo>::value)
    {
        assert("Complex(RootTwo(1, 0), RootTwo(2, 0))" == c1.toString());
        assert("Complex(RootTwo(4, 0), RootTwo(9, 0))" == c2.toString());
        std::cout << "\ttoString tests passed" << std::endl;
    }
    else
    {
        templateError(typeid(T).name(), "Complex");
    }

    assert(Complex<T>(5, 2) == c1 + 4);
    std::cout << "\tscalar sum test passed" << std::endl;

    assert(Complex<T>(-1, 2) == c1 - 2);
    std::cout << "\tscalar difference test passed" << std::endl;

    assert(Complex<T>(12, 27) == c2 * 3);
    std::cout << "\tscalar product test passed" << std::endl;

    assert(Complex<T>(5, 11) == c1 + c2);
    std::cout << "\tsum test passed" << std::endl;

    assert(Complex<T>(-3, -7) == c1 - c2);
    std::cout << "\tdifference test passed" << std::endl;

    assert(Complex<T>(-53, 2) == c2 * c3);
    std::cout << "\tproduct test passed" << std::endl;

    assert(Complex<T>(-1, -2) == -c1);
    std::cout << "\tnegation test passed" << std::endl;

    // abs doesn't do anything for complex numbers.
    assert(Complex<T>(1, 2) == c1.abs());
    assert(Complex<T>(-2, 5) == c3.abs());
    std::cout << "\tabs tests passed" << std::endl;

    // signum is always 1 for complex numbers.
    assert(1 == c1.signum());
    assert(1 == c2.signum());
    std::cout << "\tsignum tests passed" << std::endl;

    assert(Complex<T>(1, -2) == c1.adj());
    assert(Complex<T>(-2, -5) == c3.adj());
    std::cout << "\tadj tests passed" << std::endl;
}

template <typename T>
void testOmega()
{
    std::cout << "Omega<" << typeid(T).name() << "> testing:" << std::endl;

    Omega<T> o0 = Omega<T>(0, 0, 0, 0);
    Omega<T> o1 = Omega<T>(1, -3, 3, -7);
    Omega<T> o1copy = Omega<T>(1, -3, 3, -7);
    Omega<T> o2 = Omega<T>(3, 7, 12, 2);
    Omega<T> o2copy = Omega<T>(3, 7, 12, 2);
    Omega<T> o3 = Omega<T>(5, 9, -4, -9);

    assert(o1 == o1copy);
    assert(o2 == o2copy);
    assert(o1 != o2);
    assert(o2 != o3);
    assert(o1 != o3);
    assert(!(o1 != o1copy));
    assert(!(o2 != o2copy));
    assert(!(o1 == o2));
    assert(!(o2 == o3));
    assert(!(o1 == o3));
    std::cout << "\tequality tests passed" << std::endl;

    if constexpr (std::is_same<T, int>::value || std::is_same<T, Integer>::value || std::is_same<T, Rational>::value)
    {
        assert("Omega(0, 0, 0, 0)" == o0.toString());
        assert("Omega(1, -3, 3, -7)" == o1.toString());
        assert("Omega(3, 7, 12, 2)" == o2.toString());
        assert("Omega(5, 9, -4, -9)" == o3.toString());
        std::cout << "\ttoString tests passed" << std::endl;
    }
    else if constexpr (std::is_same<T, double>::value)
    {
        assert("Omega(0.000000, 0.000000, 0.000000, 0.000000)" == o0.toString());
        assert("Omega(1.000000, -3.000000, 3.000000, -7.000000)" == o1.toString());
        assert("Omega(3.000000, 7.000000, 12.000000, 2.000000)" == o2.toString());
        assert("Omega(5.000000, 9.000000, -4.000000, -9.000000)" == o3.toString());
        std::cout << "\ttoString tests passed" << std::endl;
    }
    else if constexpr (std::is_same<T, ZDyadic>::value)
    {
        assert("Omega(Dyadic(0, 0), Dyadic(0, 0), Dyadic(0, 0), Dyadic(0, 0))" == o0.toString());
        assert("Omega(Dyadic(1, 0), Dyadic(-3, 0), Dyadic(3, 0), Dyadic(-7, 0))" == o1.toString());
        assert("Omega(Dyadic(3, 0), Dyadic(7, 0), Dyadic(12, 0), Dyadic(2, 0))" == o2.toString());
        assert("Omega(Dyadic(5, 0), Dyadic(9, 0), Dyadic(-4, 0), Dyadic(-9, 0))" == o3.toString());
        std::cout << "\ttoString tests passed" << std::endl;
    }
    else
    {
        templateError(typeid(T).name(), "Omega");
    }

    assert(Omega<T>(0, 0, 0, 0) == o0.copy());
    assert(Omega<T>(1, -3, 3, -7) == o1.copy());
    std::cout << "\tcopy tests passed" << std::endl;

    assert(Omega<T>(4, 4, 15, -5) == o1 + o2);
    assert(Omega<T>(8, 16, 8, -7) == o2 + o3);
    std::cout << "\tsum tests passed" << std::endl;

    assert(Omega<T>(-2, -10, -9, -9) == o1 - o2);
    assert(Omega<T>(-2, -2, 16, 11) == o2 - o3);
    std::cout << "\tdifference tests passed" << std::endl;

    assert(Omega<T>(0, 0, 0, 0) == o0 * o2);
    assert(Omega<T>(-34, -22, -76, -14) == o1 * o2);
    assert(Omega<T>(63, -108, -178, -129) == o2 * o3);
    std::cout << "\tproduct tests passed" << std::endl;

    assert(Omega<T>(0, 0, 0, 0) == -o0);
    assert(Omega<T>(-1, 3, -3, 7) == -o1);
    assert(Omega<T>(-3, -7, -12, -2) == -o2);
    std::cout << "\tnegation tests passed" << std::endl;

    assert(Omega<T>(0, 0, 0, 0) == o0.abs());
    assert(Omega<T>(1, -3, 3, -7) == o1.abs());
    assert(Omega<T>(3, 7, 12, 2) == o2.abs());
    std::cout << "\tabs tests passed" << std::endl;

    assert(1 == o0.signum());
    assert(1 == o1.signum());
    assert(1 == o2.signum());
    assert(1 == o3.signum());
    std::cout << "\tsignum tests passed" << std::endl;

    assert(Omega<T>(0, 0, 0, 0) == o0.adj());
    assert(Omega<T>(-3, 3, -1, -7) == o1.adj());
    assert(Omega<T>(-12, -7, -3, 2) == o2.adj());
    assert(Omega<T>(4, -9, -5, -9) == o3.adj());
    std::cout << "\tadj tests passed" << std::endl;

    if constexpr (std::is_same<T, double>::value)
    {
        assert(Omega<T>(1.7114914425427872e-2, 3.0562347188264057e-2, -5.256723716381418e-2, -0.1295843520782396) == o1.recip());
        assert(Omega<T>(-0.11208737066841847, 3.3092461816390216e-2, -2.4634586960091973e-4, -5.7070126457546395e-2) == o2.recip());
        assert(Omega<T>(3.2468311407893156e-2, -6.945499620136751e-2, 2.419129113519133e-2, -4.066536047023071e-2) == o3.recip());
        std::cout << "\trecip tests passed" << std::endl;
    }
    else if constexpr (std::is_same<T, Rational>::value)
    {
        assert(Omega<T>(7_mpq / 409, 25_mpq / 818, -43_mpq / 818, -53_mpq / 409) == o1.recip());
        assert(Omega<T>(-1365_mpq / 12178, 403_mpq / 12178, -3_mpq / 12178, -695_mpq / 12178) == o2.recip());
        assert(Omega<T>(812_mpq / 25009, -1737_mpq / 25009, 605_mpq / 25009, -1017_mpq / 25009) == o3.recip());
        std::cout << "\trecip tests passed" << std::endl;
    }
    else if constexpr (std::is_same<T, int>::value || std::is_same<T, Integer>::value || std::is_same<T, ZDyadic>::value)
    {
        // recip is not defined when T = int, T = Integer, or T = ZDyadic.
    }
    else
    {
        templateError(typeid(T).name(), "Omega");
    }
}

int main()
{
    testUtilityFunctions();
    std::cout << std::endl;
    testTypeConversions();
    std::cout << std::endl;

    testHalfRing();
    std::cout << std::endl;
    testRootTwoRing();
    std::cout << std::endl;
    testRootHalfRing();
    std::cout << std::endl;
    testComplexRing();
    std::cout << std::endl;
    testOmegaRing();
    std::cout << std::endl;
    testAdjoint();
    std::cout << std::endl;
    testAdjoint2();
    std::cout << std::endl;
    testNormedRing();
    std::cout << std::endl;
    testFloor();
    std::cout << std::endl;

    testRootTwoIntegral<int>();
    std::cout << std::endl;
    testRootTwoIntegral<Integer>();
    std::cout << std::endl;
    testRootTwoDyadic<int>();
    std::cout << std::endl;
    testRootTwoDyadic<Integer>();
    std::cout << std::endl;
    testRootTwoRational();
    std::cout << std::endl;

    testComplex<int>();
    std::cout << std::endl;
    testComplex<Integer>();
    std::cout << std::endl;
    testComplex<double>();
    std::cout << std::endl;
    testComplex<Rational>();
    std::cout << std::endl;
    testComplex<ZDyadic>();
    std::cout << std::endl;
    testComplex<DRootTwo>();
    std::cout << std::endl;
    testComplex<QRootTwo>();
    std::cout << std::endl;

    testDyadic<int>();
    std::cout << std::endl;
    testDyadic<Integer>();
    std::cout << std::endl;

    testZ2();
    std::cout << std::endl;

    testOmega<int>();
    std::cout << std::endl;
    testOmega<Integer>();
    std::cout << std::endl;
    testOmega<double>();
    std::cout << std::endl;
    testOmega<Rational>();
    std::cout << std::endl;
    testOmega<ZDyadic>();
    std::cout << std::endl;

    return 0;
}
