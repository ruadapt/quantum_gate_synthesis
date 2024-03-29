#define BOOST_TEST_MODULE ring
#include "comparisons.h"
#include "../src/ring.h"
#include <boost/test/included/unit_test.hpp>
#include <assert.h>
#include <iostream>

void templateError(std::string T, std::string baseType)
{
    std::cout << "Unexpected template type " << T << " for " << baseType << std::endl;
    exit(1);
}

BOOST_AUTO_TEST_CASE(test_utility_functions)
{
    assert(1 == ring::shift<int>(1, 0));
    assert(1 == ring::shift<Integer>(1, 0));
    assert(5120 == ring::shift<int>(5, 10));
    assert(5120 == ring::shift<Integer>(5, 10));
    assert(8034690221294951377709810461705813012611014968913964176506880_mpz == ring::shift<Integer>(5, 200));
    assert(540412499179089513452888265605116091644885096989720576_mpz ==
           ring::shift<Integer>(515377520732011331036461129765621272702107522001_mpz, 20));
    // Fails at runtime because the shift ammount is too large.
    // assert(1 == ring::shift<Integer>(3, 515377520732011331036461129765621272702107522001_mpz));

    // shiftL should be the same as shift.
    assert(1 == ring::shiftL<int>(1, 0));
    assert(1 == ring::shiftL<Integer>(1, 0));
    assert(5120 == ring::shiftL<int>(5, 10));
    assert(5120 == ring::shiftL<Integer>(5, 10));
    assert(8034690221294951377709810461705813012611014968913964176506880_mpz == ring::shiftL<Integer>(5, 200));
    assert(540412499179089513452888265605116091644885096989720576_mpz ==
           ring::shiftL<Integer>(515377520732011331036461129765621272702107522001_mpz, 20));

    assert(1 == ring::shiftR<int>(1, 0));
    assert(1 == ring::shiftR<Integer>(1, 0));
    assert(5 == ring::shiftR<int>(5120, 10));
    assert(5 == ring::shiftR<Integer>(5120, 10));
    assert(5 == ring::shiftR<Integer>(8034690221294951377709810461705813012611014968913964176506880_mpz, 200));
    assert(515377520732011331036461129765621272702107522001_mpz ==
           ring::shiftR<Integer>(540412499179089513452888265605116091644885096989720576_mpz, 20));

    assert(1 == ring::exp2<int>(0));
    assert(1 == ring::exp2<Integer>(0));
    assert(1024 == ring::exp2<int>(10));
    assert(1024 == ring::exp2<Integer>(10));
    assert(1606938044258990275541962092341162602522202993782792835301376_mpz == ring::exp2<Integer>(200));

    assert(2 == ring::intsqrt(4));
    assert(11096085937082_mpz == ring::intsqrt(123123123123123123123123123_mpz));
    assert(456456456456456456_mpz == ring::intsqrt(208352496640784928656512368224079936_mpz));

    assert(!ring::log2(0).has_value());
    assert(0 == ring::log2(1).value());
    assert(!ring::log2(755).has_value());
    assert(10 == ring::log2(1024).value());
    assert(200 == ring::log2(1606938044258990275541962092341162602522202993782792835301376_mpz).value());
    assert(!ring::log2(1606938044258990275541962092341162602522202993782792835301377_mpz).has_value());

    assert(QRootTwo(32) == ring::pow_non_neg(QRootTwo(2), 5));
    assert(ZOmega(0) == ring::pow_non_neg(ZOmega(0), 5));
    assert(ZDyadic(1099511627776_mpz) == ring::pow_non_neg(ZDyadic(2), 40));
    assert(1 == ring::pow_non_neg(3, 0));
    assert(1000_mpz == ring::pow_non_neg<Integer>(10_mpz, 3));

    assert(QRootTwo(32) == ring::pow_int(QRootTwo(2), 5));
    assert(QRootTwo(1_mpq / 27, 0) == ring::pow_int(QRootTwo(3), -3));

    assert(3_mpz == ring::div(10_mpz, 3_mpz));
    assert(3_mpz == ring::div(-10_mpz, -3_mpz));
    assert(-4_mpz == ring::div(10_mpz, -3_mpz));
    assert(-4_mpz == ring::div(-10_mpz, 3_mpz));
    assert(3 == ring::div(10, 3));
    assert(3 == ring::div(-10, -3));
    assert(-4 == ring::div(10, -3));
    assert(-4 == ring::div(-10, 3));
}

BOOST_AUTO_TEST_CASE(test_type_conversions)
{
    Real d = Real(2);
    assert(2.0 == d);
    Integer i = Integer(2);
    assert(2_mpz == i);
    Rational r = Rational(2);
    assert(2_mpq == r);
    Dyadic<int> dInt = Dyadic<int>(2);
    assert(Dyadic<int>(2, 0) == dInt);
    ZDyadic dInteger = ZDyadic(2);
    assert(ZDyadic(2_mpz, 0_mpz) == dInteger);
    RootTwo<int> rInt = RootTwo<int>(2);
    assert(RootTwo<int>(2, 0) == rInt);
    ZRootTwo rInteger = ZRootTwo(2);
    assert(ZRootTwo(2_mpz, 0_mpz) == rInteger);
    RootTwo<Real> rReal = RootTwo<Real>(2);
    assert(RootTwo<Real>(2.0, 0.0) == rReal);
    QRootTwo rRational = QRootTwo(2);
    assert(QRootTwo(2_mpq, 0_mpq) == rRational);
    RootTwo<Dyadic<int>> rDyadicInt = RootTwo<Dyadic<int>>(2);
    assert(RootTwo<Dyadic<int>>(Dyadic<int>(2, 0), Dyadic<int>(0, 0)) == rDyadicInt);
    DRootTwo rDyadicInteger = DRootTwo(2);
    assert(DRootTwo(ZDyadic(2_mpz, 0_mpz), ZDyadic(0_mpz, 0_mpz)) == rDyadicInteger);
    QRComplex cQRootTwo = QRComplex(2);
    assert(QRComplex(QRootTwo(2_mpq, 0_mpq), QRootTwo(0_mpq, 0_mpq)) == cQRootTwo);
    DRComplex cDRootTwo = DRComplex(2);
    assert(DRComplex(DRootTwo(ZDyadic(2, 0), ZDyadic(0, 0)), DRootTwo(ZDyadic(0, 0), ZDyadic(0, 0))) == cDRootTwo);
    Omega<int> oInt = Omega<int>(2);
    assert(Omega<int>(0, 0, 0, 2) == oInt);
    ZOmega oInteger = ZOmega(2);
    assert(ZOmega(0_mpz, 0_mpz, 0_mpz, 2_mpz) == oInteger);
    QOmega oRational = QOmega(2);
    assert(QOmega(0_mpq, 0_mpq, 0_mpq, 2_mpq) == oRational);
    DOmega oDyadicInteger = DOmega(2);
    assert(DOmega(ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(2, 0)) == oDyadicInteger);

    Integer n = 123412341234123412341234123412341234123412341234123412341234123412341234_mpz;

    // Make a new variable scope to avoid name collisions.
    {
        Integer i = Integer(n);
        assert(n == i);
        Rational r = Rational(n);
        assert(n == r);
        ZDyadic dInteger = ZDyadic(n);
        assert(ZDyadic(n, 0) == dInteger);
        ZRootTwo rInteger = ZRootTwo(n);
        assert(ZRootTwo(n, 0) == rInteger);
        QRootTwo rRational = QRootTwo(n);
        assert(QRootTwo(n, 0) == rRational);
        DRootTwo rDyadicInteger = DRootTwo(n);
        assert(DRootTwo(ZDyadic(n, 0), ZDyadic(0, 0)) == rDyadicInteger);
        QRComplex cQRootTwo = QRComplex(n);
        assert(QRComplex(QRootTwo(n, 0), QRootTwo(0, 0)) == cQRootTwo);
        DRComplex cDRootTwo = DRComplex(n);
        assert(DRComplex(DRootTwo(ZDyadic(n, 0), ZDyadic(0, 0)), DRootTwo(ZDyadic(0, 0), ZDyadic(0, 0))) == cDRootTwo);
        ZOmega oInteger = ZOmega(n);
        assert(ZOmega(0, 0, 0, n) == oInteger);
        QOmega oRational = QOmega(n);
        assert(QOmega(0, 0, 0, n) == oRational);
        DOmega oDyadicInteger = DOmega(n);
        assert(DOmega(ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(n, 0)) == oDyadicInteger);
    }
    return;

    assert(QRootTwo(2) == ring::fromQRootTwo<QRootTwo>(2));
    assert(QRComplex(2) == ring::fromQRootTwo<QRComplex>(2));
    assert(QOmega(2) == ring::fromQRootTwo<QOmega>(2));

    assert(ZComplex(2) == ring::fromZComplex<ZComplex>(2));
    assert(DComplex(2) == ring::fromZComplex<DComplex>(2));
    assert(QComplex(2) == ring::fromZComplex<QComplex>(2));
    assert(DRComplex(2) == ring::fromZComplex<DRComplex>(2));
    assert(QRComplex(2) == ring::fromZComplex<QRComplex>(2));
    assert(CReal(2) == ring::fromZComplex<CReal>(2));
    assert(ZOmega(2) == ring::fromZComplex<ZOmega>(2));
    assert(DOmega(2) == ring::fromZComplex<DOmega>(2));
    assert(QOmega(2) == ring::fromZComplex<QOmega>(2));

    assert(DComplex(2) == ring::fromDComplex<DComplex>(2));
    assert(QComplex(2) == ring::fromDComplex<QComplex>(2));
    assert(DRComplex(2) == ring::fromDComplex<DRComplex>(2));
    assert(QRComplex(2) == ring::fromDComplex<QRComplex>(2));
    assert(CReal(2) == ring::fromDComplex<CReal>(2));
    assert(DOmega(2) == ring::fromDComplex<DOmega>(2));
    assert(QOmega(2) == ring::fromDComplex<QOmega>(2));

    assert(QComplex(2) == ring::fromQComplex<QComplex>(2));
    assert(QRComplex(2) == ring::fromQComplex<QRComplex>(2));
    assert(CReal(2) == ring::fromQComplex<CReal>(2));
    assert(QOmega(2) == ring::fromQComplex<QOmega>(2));

    assert(DRComplex(2) == ring::fromDRComplex<DRComplex>(2));
    assert(QRComplex(2) == ring::fromDRComplex<QRComplex>(2));
    assert(CReal(2) == ring::fromDRComplex<CReal>(2));
    assert(DOmega(2) == ring::fromDRComplex<DOmega>(2));
    assert(QOmega(2) == ring::fromDRComplex<QOmega>(2));

    assert(QRComplex(2) == ring::fromQRComplex<QRComplex>(2));
    assert(CReal(2) == ring::fromQRComplex<CReal>(2));
    assert(QOmega(2) == ring::fromQRComplex<QOmega>(2));

    assert(DRComplex(2) == ring::fromZOmega<DRComplex>(2));
    assert(QRComplex(2) == ring::fromZOmega<QRComplex>(2));
    assert(CReal(2) == ring::fromZOmega<CReal>(2));
    assert(ZOmega(2) == ring::fromZOmega<ZOmega>(2));
    assert(DOmega(2) == ring::fromZOmega<DOmega>(2));
    assert(QOmega(2) == ring::fromZOmega<QOmega>(2));

    assert(DRComplex(2) == ring::fromDOmega<DRComplex>(2));
    assert(QRComplex(2) == ring::fromDOmega<QRComplex>(2));
    assert(CReal(2) == ring::fromDOmega<CReal>(2));
    assert(DOmega(2) == ring::fromDOmega<DOmega>(2));
    assert(QOmega(2) == ring::fromDOmega<QOmega>(2));

    assert(QRComplex(2) == ring::fromQOmega<QRComplex>(2));
    assert(CReal(2) == ring::fromQOmega<CReal>(2));
    assert(QOmega(2) == ring::fromQOmega<QOmega>(2));
}

BOOST_AUTO_TEST_CASE(test_fractional)
{
    RootTwo<Real> r = RootTwo<Real>(1.25, -2.45).recip();
    assert(approx_equal(-0.119703136222169, r.a()));
    assert(approx_equal(-0.23461814699545125, r.b()));
    assert(QRootTwo(-431061208_mpq / 951316543, 923605504_mpq / 951316543) ==
           QRootTwo(123_mpq / 456, 456_mpq / 789).recip());

    assert(QComplex(431061208_mpq / 650067905, -923605504_mpq / 650067905) ==
           QComplex(123_mpq / 456, 456_mpq / 789).recip());

    Omega<Real> o1 = Omega<Real>(1, -3, 3, -7).recip();
    assert(approx_equal(1.7114914425427872e-2, o1.a()));
    assert(approx_equal(3.0562347188264057e-2, o1.b()));
    assert(approx_equal(-5.256723716381418e-2, o1.c()));
    assert(approx_equal(-0.1295843520782396, o1.d()));

    Omega<Real> o2 = Omega<Real>(3, 7, 12, 2).recip();
    assert(approx_equal(-0.11208737066841847, o2.a()));
    assert(approx_equal(3.3092461816390216e-2, o2.b()));
    assert(approx_equal(-2.4634586960091973e-4, o2.c()));
    assert(approx_equal(-5.7070126457546395e-2, o2.d()));

    Omega<Real> o3 = Omega<Real>(5, 9, -4, -9).recip();
    assert(approx_equal(3.2468311407893156e-2, o3.a()));
    assert(approx_equal(-6.945499620136751e-2, o3.b()));
    assert(approx_equal(2.419129113519133e-2, o3.c()));
    assert(approx_equal(-4.066536047023071e-2, o3.d()));

    assert(QOmega(7_mpq / 409, 25_mpq / 818, -43_mpq / 818, -53_mpq / 409) ==
           QOmega(1, -3, 3, -7).recip());
    assert(QOmega(-1365_mpq / 12178, 403_mpq / 12178, -3_mpq / 12178, -695_mpq / 12178) ==
           QOmega(3, 7, 12, 2).recip());
    assert(QOmega(812_mpq / 25009, -1737_mpq / 25009, 605_mpq / 25009, -1017_mpq / 25009) ==
           QOmega(5, 9, -4, -9).recip());

    assert(QRootTwo(123_mpq / 456, 0) == ring::fromRational<QRootTwo>(123_mpq / 456));
    assert(QComplex(123_mpq / 456, 0) == ring::fromRational<QComplex>(123_mpq / 456));
    assert(QOmega(0, 0, 0, 41_mpq / 152) == ring::fromRational<QOmega>(123_mpq / 456));
}

BOOST_AUTO_TEST_CASE(test_half_ring)
{
    assert(0.5 == ring::half<Real>());
    assert(1_mpq / 2 == ring::half<Rational>());
    assert(Dyadic<int>(1, 1) == ring::half<Dyadic<int>>());
    assert(ZDyadic(1, 1) == ring::half<ZDyadic>());
    assert(QRootTwo(1_mpq / 2, 0_mpq) == ring::half<QRootTwo>());
    assert(CReal(0.5, 0) == ring::half<CReal>());
    assert(Omega<Real>(0, 0, 0, 0.5) == ring::half<Omega<Real>>());
    assert(QOmega(0_mpq, 0_mpq, 0_mpq, 1_mpq / 2) == ring::half<QOmega>());
    assert(DOmega(ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(1, 1)) == ring::half<DOmega>());

    assert(0.625 == ring::fromDyadic<Real>(Dyadic<int>(5, 3)));
    assert(40 == ring::fromDyadic<Real>(Dyadic<int>(5, -3)));
    assert(5_mpq / 8 == ring::fromDyadic<Rational>(Dyadic<int>(5, 3)));
    assert(40_mpq == ring::fromDyadic<Rational>(Dyadic<int>(5, -3)));
    assert(QRootTwo(5_mpq / 8, 0_mpq) == ring::fromDyadic<QRootTwo>(Dyadic<int>(5, 3)));
    assert(QRootTwo(40_mpq, 0_mpq) == ring::fromDyadic<QRootTwo>(Dyadic<int>(5, -3)));
    assert(CReal(0.625, 0) == ring::fromDyadic<CReal>(Dyadic<int>(5, 3)));
    assert(CReal(40, 0) == ring::fromDyadic<CReal>(Dyadic<int>(5, -3)));
    assert(Omega<Real>(0, 0, 0, 0.625) == ring::fromDyadic<Omega<Real>>(Dyadic<int>(5, 3)));
    assert(Omega<Real>(0, 0, 0, 40) == ring::fromDyadic<Omega<Real>>(Dyadic<int>(5, -3)));
    assert(QOmega(0_mpq, 0_mpq, 0_mpq, 5_mpq / 8) ==
           ring::fromDyadic<QOmega>(Dyadic<int>(5, 3)));
    assert(QOmega(0_mpq, 0_mpq, 0_mpq, 40_mpq) ==
           ring::fromDyadic<QOmega>(Dyadic<int>(5, -3)));
    assert(DOmega(ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(5, 3)) ==
           ring::fromDyadic<DOmega>(Dyadic<int>(5, 3)));
    assert(DOmega(ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(5, -3)) ==
           ring::fromDyadic<DOmega>(Dyadic<int>(5, -3)));

    assert(0.625 == ring::fromDyadic<Real>(ZDyadic(5, 3)));
    assert(40 == ring::fromDyadic<Real>(ZDyadic(5, -3)));
    assert(5_mpq / 8 == ring::fromDyadic<Rational>(ZDyadic(5, 3)));
    assert(40_mpq == ring::fromDyadic<Rational>(ZDyadic(5, -3)));
    assert(QRootTwo(5_mpq / 8, 0_mpq) == ring::fromDyadic<QRootTwo>(ZDyadic(5, 3)));
    assert(QRootTwo(40_mpq, 0_mpq) == ring::fromDyadic<QRootTwo>(ZDyadic(5, -3)));
    assert(CReal(0.625, 0) == ring::fromDyadic<CReal>(ZDyadic(5, 3)));
    assert(CReal(40, 0) == ring::fromDyadic<CReal>(ZDyadic(5, -3)));
    assert(Omega<Real>(0, 0, 0, 0.625) == ring::fromDyadic<Omega<Real>>(ZDyadic(5, 3)));
    assert(Omega<Real>(0, 0, 0, 40) == ring::fromDyadic<Omega<Real>>(ZDyadic(5, -3)));
    assert(QOmega(0_mpq, 0_mpq, 0_mpq, 5_mpq / 8) ==
           ring::fromDyadic<QOmega>(ZDyadic(5, 3)));
    assert(QOmega(0_mpq, 0_mpq, 0_mpq, 40_mpq) ==
           ring::fromDyadic<QOmega>(ZDyadic(5, -3)));
    assert(DOmega(ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(5, 3)) ==
           ring::fromDyadic<DOmega>(ZDyadic(5, 3)));
    assert(DOmega(ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(0, 0), ZDyadic(5, -3)) ==
           ring::fromDyadic<DOmega>(ZDyadic(5, -3)));
}

BOOST_AUTO_TEST_CASE(test_roottwo_ring)
{
    assert(approx_equal(1.4142135623730951, ring::roottwo<Real>()));
    assert(approx_equal(1.4142135623730951, ring::roottwo<CReal>().a()));
    assert(0 == ring::roottwo<CReal>().b());
    assert(ZRootTwo(0, 1) == ring::roottwo<ZRootTwo>());
    assert(QRootTwo(0_mpq, 1_mpq) == ring::roottwo<QRootTwo>());
    assert(ZOmega(-1, 0, 1, 0) == ring::roottwo<ZOmega>());
    assert(DOmega(-1, 0, 1, 0) == ring::roottwo<DOmega>());
    assert(QOmega(-1, 0, 1, 0) == ring::roottwo<QOmega>());

    ZRootTwo z = ZRootTwo(10, -7);
    assert(approx_equal(0.10050506338833465, ring::fromZRootTwo<Real>(z)));
    assert(approx_equal(0.10050506338833465, ring::fromZRootTwo<CReal>(z).a()));
    assert(0 == ring::fromZRootTwo<CReal>(z).b());
    assert(ZRootTwo(10, -7) == ring::fromZRootTwo<ZRootTwo>(z));
    assert(QRootTwo(10_mpq, -7_mpq) == ring::fromZRootTwo<QRootTwo>(z));
    assert(ZOmega(7, 0, -7, 10) == ring::fromZRootTwo<ZOmega>(z));
    assert(DOmega(7, 0, -7, 10) == ring::fromZRootTwo<DOmega>(z));
    assert(QOmega(7, 0, -7, 10) == ring::fromZRootTwo<QOmega>(z));
}

BOOST_AUTO_TEST_CASE(test_roothalf_ring)
{
    assert(approx_equal(0.7071067811865476, ring::roothalf<Real>()));
    assert(approx_equal(0.7071067811865476, ring::roothalf<CReal>().a()));
    assert(0 == ring::roothalf<CReal>().b());
    assert(QRootTwo(0, 1_mpq / 2) == ring::roothalf<QRootTwo>());
    assert(DOmega(ZDyadic(-1, 1), 0, ZDyadic(1, 1), 0) == ring::roothalf<DOmega>());
    assert(QOmega(-1_mpq / 2, 0, 1_mpq / 2, 0) == ring::roothalf<QOmega>());

    DRootTwo d = DRootTwo(ZDyadic(5, 2), ZDyadic(-3, 5));
    assert(approx_equal(1.1174174785275224, ring::fromDRootTwo<Real>(d)));
    assert(approx_equal(1.1174174785275224, ring::fromDRootTwo<CReal>(d).a()));
    assert(0 == ring::fromDRootTwo<CReal>(d).b());
    assert(QRootTwo(5_mpq / 4, -3_mpq / 32) == ring::fromDRootTwo<QRootTwo>(d));
    assert(DOmega(ZDyadic(3, 5), 0, ZDyadic(-3, 5), ZDyadic(5, 2)) == ring::fromDRootTwo<DOmega>(d));
    assert(QOmega(3_mpq / 32, 0, -3_mpq / 32, 5_mpq / 4) == ring::fromDRootTwo<QOmega>(d));
}

BOOST_AUTO_TEST_CASE(test_complex_ring)
{

    assert(Complex<int>(0, 1) == ring::i<Complex<int>>());
    assert(ZComplex(0, 1) == ring::i<ZComplex>());
    assert(Complex<Complex<int>>(Complex<int>(0, 0), Complex<int>(1, 0)) == ring::i<Complex<Complex<int>>>());
    assert(RootTwo<Complex<int>>(Complex<int>(0, 1), Complex<int>(0, 0)) == ring::i<RootTwo<Complex<int>>>());
    assert(RootTwo<ZComplex>(ZComplex(0, 1), ZComplex(0, 0)) == ring::i<RootTwo<ZComplex>>());
    assert(ZOmega(0, 1, 0, 0) == ring::i<ZOmega>());
    assert(DOmega(0, 1, 0, 0) == ring::i<DOmega>());
    assert(QOmega(0, 1, 0, 0) == ring::i<QOmega>());
}

BOOST_AUTO_TEST_CASE(test_omega_ring)
{

    assert(approx_equal(0.7071067811865476, ring::omega<CReal>().a()));
    assert(approx_equal(0.7071067811865476, ring::omega<CReal>().b()));
    assert(DRComplex(DRootTwo(0, ZDyadic(1, 1)), DRootTwo(0, ZDyadic(1, 1))) ==
           ring::omega<DRComplex>());
    assert(QRComplex(QRootTwo(0, 1_mpq / 2), QRootTwo(0, 1_mpq / 2)) ==
           ring::omega<QRComplex>());

    assert(RootTwo<DComplex>(DComplex(0, 0), DComplex(ZDyadic(1, 1), ZDyadic(1, 1))) ==
           ring::omega<RootTwo<DComplex>>());

    assert(ZOmega(0, 0, 1, 0) == ring::omega<ZOmega>());
    assert(DOmega(0, 0, 1, 0) == ring::omega<DOmega>());
    assert(QOmega(0, 0, 1, 0) == ring::omega<QOmega>());
}

BOOST_AUTO_TEST_CASE(test_normed_ring)
{

    assert(10_mpz == ring::norm<int>(10));
    assert(10_mpz == ring::norm<Integer>(10));
    assert(85_mpz == ring::norm<Complex<int>>(Complex<int>(6, -7)));
    assert(85_mpz == ring::norm<ZComplex>(ZComplex(6, -7)));
    assert(-62_mpz == ring::norm<RootTwo<int>>(RootTwo<int>(6, -7)));
    assert(-62_mpz == ring::norm<ZRootTwo>(ZRootTwo(6, -7)));
    // norm = (3^2 + 4^2)^2 - 2 * (5^2 + 2^2)^2 = 35^2 - 2 * 29^2 = -1057.
    assert(-1057_mpz == ring::norm<RootTwo<Complex<int>>>(
                            RootTwo<Complex<int>>(Complex<int>(3, 4), Complex<int>(5, 2))));
    // norm = (3^2 - 2 * 4^2)^2 + (5^2 - 2 * 2^2)^2 = (-23)^2 + (17)^2 = 818.
    assert(818_mpq == ring::norm<Complex<RootTwo<int>>>(
                          Complex<RootTwo<int>>(RootTwo<int>(3, 4), RootTwo<int>(5, 2))));
    assert(6562_mpq == ring::norm<Omega<int>>(Omega<int>(7, 3, -2, 6)));
    assert(6562_mpq == ring::norm<ZOmega>(ZOmega(7, 3, -2, 6)));
}

BOOST_AUTO_TEST_CASE(test_adjoint)
{

    assert(7 == ring::adj<int>(7));
    assert(123123123123123123_mpz == ring::adj<Integer>(123123123123123123_mpz));
    assert(0.123 == ring::adj<Real>(0.123));
    assert(123_mpq / 456 == ring::adj<Rational>(123_mpq / 456));
    assert(Z2(0) == ring::adj<Z2>(0));
    assert(Z2(1) == ring::adj<Z2>(1));
    assert(Dyadic<int>(2, 7) == ring::adj<Dyadic<int>>(Dyadic<int>(2, 7)));
    assert(ZDyadic(2, 7) == ring::adj<ZDyadic>(ZDyadic(2, 7)));
    assert(Complex<int>(3, -4) == ring::adj<Complex<int>>(Complex<int>(3, 4)));
    assert(ZComplex(3, -4) == ring::adj<ZComplex>(ZComplex(3, 4)));
    assert(RootTwo<int>(4, 7) == ring::adj<RootTwo<int>>(RootTwo<int>(4, 7)));
    assert(ZRootTwo(4, 7) == ring::adj<ZRootTwo>(ZRootTwo(4, 7)));
    assert(
        RootTwo<Complex<int>>(Complex<int>(4, -5), Complex<int>(6, -7)) ==
        ring::adj<RootTwo<Complex<int>>>(RootTwo<Complex<int>>(Complex<int>(4, 5), Complex<int>(6, 7))));
    assert(
        RootTwo<ZComplex>(ZComplex(4, -5), ZComplex(6, -7)) ==
        ring::adj<RootTwo<ZComplex>>(RootTwo<ZComplex>(ZComplex(4, 5), ZComplex(6, 7))));
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
}

BOOST_AUTO_TEST_CASE(test_adjoint2)
{

    assert(7 == ring::adj2<int>(7));
    assert(123123123123123123_mpz == ring::adj2<Integer>(123123123123123123_mpz));
    assert(123_mpq / 456 == ring::adj2<Rational>(123_mpq / 456));
    assert(Z2(0) == ring::adj2<Z2>(0));
    assert(Z2(1) == ring::adj2<Z2>(1));
    assert(Dyadic<int>(2, 7) == ring::adj2<Dyadic<int>>(Dyadic<int>(2, 7)));
    assert(ZDyadic(2, 7) == ring::adj2<ZDyadic>(ZDyadic(2, 7)));
    assert(Complex<int>(3, 4) == ring::adj2<Complex<int>>(Complex<int>(3, 4)));
    assert(ZComplex(3, 4) == ring::adj2<ZComplex>(ZComplex(3, 4)));
    assert(RootTwo<int>(4, -7) == ring::adj2<RootTwo<int>>(RootTwo<int>(4, 7)));
    assert(ZRootTwo(4, -7) == ring::adj2<ZRootTwo>(ZRootTwo(4, 7)));
    assert(
        RootTwo<Complex<int>>(Complex<int>(4, 5), Complex<int>(-6, -7)) ==
        ring::adj2<RootTwo<Complex<int>>>(RootTwo<Complex<int>>(Complex<int>(4, 5), Complex<int>(6, 7))));
    assert(
        RootTwo<ZComplex>(ZComplex(4, 5), ZComplex(-6, -7)) ==
        ring::adj2<RootTwo<ZComplex>>(RootTwo<ZComplex>(ZComplex(4, 5), ZComplex(6, 7))));
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
}

BOOST_AUTO_TEST_CASE(test_floor)
{

    // Real
    assert(5 == ring::floor_of(Real(5.2)));
    assert(4 == ring::floor_of(Real(4)));
    assert(-3 == ring::floor_of(Real(-2.5)));
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

    // Real
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
}

BOOST_AUTO_TEST_CASE(test_to_dyadic)
{

    assert(ZDyadic(5, 8) == (ring::maybe_dyadic<ZDyadic, ZDyadic>(ZDyadic(5, 8))));
    assert(ZDyadic(5, 8) == (ring::maybe_dyadic<Rational, ZDyadic>(5_mpq / 256)));
    assert(!(ring::maybe_dyadic<Rational, ZDyadic>(5_mpq / 252).has_value()));

    assert(DRootTwo(ZDyadic(3, 2), ZDyadic(1, 3)) ==
           (ring::maybe_dyadic<Rational, ZDyadic>(QRootTwo(3_mpq / 4, 1_mpq / 8))));
    assert(!(ring::maybe_dyadic<Rational, ZDyadic>(QRootTwo(3_mpq / 5, 1_mpq / 8)).has_value()));
    assert(!(ring::maybe_dyadic<Rational, ZDyadic>(QRootTwo(3_mpq / 4, 1_mpq / 7)).has_value()));

    assert(DComplex(ZDyadic(3, 2), ZDyadic(1, 3)) ==
           (ring::maybe_dyadic<Rational, ZDyadic>(QComplex(3_mpq / 4, 1_mpq / 8))));
    assert(!(ring::maybe_dyadic<Rational, ZDyadic>(QComplex(3_mpq / 5, 1_mpq / 8)).has_value()));
    assert(!(ring::maybe_dyadic<Rational, ZDyadic>(QComplex(3_mpq / 4, 1_mpq / 7)).has_value()));

    assert(DOmega(ZDyadic(3, 2), ZDyadic(1, 3), 2, 3) ==
           (ring::maybe_dyadic<Rational, ZDyadic>(QOmega(3_mpq / 4, 1_mpq / 8, 2, 3))));
    assert(!(ring::maybe_dyadic<Rational, ZDyadic>(QOmega(3_mpq / 5, 1_mpq / 8, 2, 2)).has_value()));
    assert(!(ring::maybe_dyadic<Rational, ZDyadic>(QOmega(3_mpq / 4, 1_mpq / 7, 2, 2)).has_value()));
    assert(!(ring::maybe_dyadic<Rational, ZDyadic>(QOmega(2, 2, 3_mpq / 7, 2)).has_value()));
    assert(!(ring::maybe_dyadic<Rational, ZDyadic>(QOmega(2, 2, 2, 3_mpq / 7)).has_value()));

    assert(ZDyadic(5, 8) == (ring::to_dyadic<ZDyadic, ZDyadic>(ZDyadic(5, 8))));
    assert(ZDyadic(5, 8) == (ring::to_dyadic<Rational, ZDyadic>(5_mpq / 256)));
    assert(DRootTwo(ZDyadic(3, 2), ZDyadic(1, 3)) ==
           (ring::to_dyadic<Rational, ZDyadic>(QRootTwo(3_mpq / 4, 1_mpq / 8))));
    assert(DComplex(ZDyadic(3, 2), ZDyadic(1, 3)) ==
           (ring::to_dyadic<Rational, ZDyadic>(QComplex(3_mpq / 4, 1_mpq / 8))));
    assert(DOmega(ZDyadic(3, 2), ZDyadic(1, 3), 2, 3) ==
           (ring::to_dyadic<Rational, ZDyadic>(QOmega(3_mpq / 4, 1_mpq / 8, 2, 3))));
}

BOOST_AUTO_TEST_CASE(test_real_part)
{

    assert(2 == ring::real(Complex<int>(2, 3)));
    assert(2 == ring::real(ZComplex(2, 3)));
    assert(ZDyadic(2) == ring::real(DComplex(2, 3)));
    assert(QRootTwo(2) == ring::real(QRComplex(2, 3)));
    assert(QRootTwo(7_mpq / 8, 1_mpq / 6) == ring::real(QOmega(1_mpq / 2, 3_mpq / 4, 5_mpq / 6, 7_mpq / 8)));
    assert(DRootTwo(ZDyadic(7, 8), ZDyadic(-11, 7)) ==
           ring::real(DOmega(ZDyadic(1, 2), ZDyadic(3, 4), ZDyadic(5, 6), ZDyadic(7, 8))));
}

BOOST_AUTO_TEST_CASE(test_whole_part)
{

    assert(ZDyadic(5) == (ring::from_whole<ZDyadic, Integer>(5)));
    assert(DOmega(1, 2, 3, 4) == (ring::from_whole<DOmega, ZOmega>(ZOmega(1, 2, 3, 4))));
    assert(DRootTwo(1, 2) == (ring::from_whole<DRootTwo, ZRootTwo>(ZRootTwo(1, 2))));

    assert(5 == (ring::to_whole<ZDyadic, Integer>(ZDyadic(5))));
    assert(ZOmega(1, 2, 3, 4) == (ring::to_whole<DOmega, ZOmega>(DOmega(1, 2, 3, 4))));
    assert(ZRootTwo(1, 2) == (ring::to_whole<DRootTwo, ZRootTwo>(DRootTwo(1, 2))));
}

BOOST_AUTO_TEST_CASE(test_denom_exp)
{

    assert(21 == ring::denomexp(DRootTwo(ZDyadic(3, 7), ZDyadic(4, 13))));
    assert(34 == ring::denomexp(DOmega(ZDyadic(12, 3), ZDyadic(15, 17), ZDyadic(13, 2), ZDyadic(5, 9))));
    assert(18 == ring::denomexp(
                     DRComplex(DRootTwo(ZDyadic(1, 9), ZDyadic(2, 8)), DRootTwo(ZDyadic(3, 7), ZDyadic(4, 6)))));

    assert(DRootTwo(ZDyadic(64, 13), ZDyadic(1536, 13)) ==
           ring::denomexp_factor(DRootTwo(ZDyadic(3, 7), ZDyadic(4, 13)), 7));
    assert(DOmega(ZDyadic(12582912, 17), ZDyadic(960, 17), ZDyadic(27262976, 17), ZDyadic(81920, 17)) ==
           ring::denomexp_factor(DOmega(ZDyadic(12, 3), ZDyadic(15, 17), ZDyadic(13, 2), ZDyadic(5, 9)), 12));
    assert(DRComplex(DRootTwo(ZDyadic(32, 9), ZDyadic(4, 9)), DRootTwo(ZDyadic(64, 7), ZDyadic(12, 7))) ==
           ring::denomexp_factor(DRComplex(DRootTwo(ZDyadic(1, 9), ZDyadic(2, 8)), DRootTwo(ZDyadic(3, 7), ZDyadic(4, 6))), 5));
}

BOOST_AUTO_TEST_CASE(test_to_qomega)
{

    assert(QOmega(5) == ring::toQOmega<Integer>(5));
    assert(QOmega(0, 0, 0, 5_mpq / 6) == ring::toQOmega<Rational>(5_mpq / 6));
    assert(QOmega(0, 0, 0, 5_mpq / 64) == ring::toQOmega<Dyadic<int>>(Dyadic<int>(5, 6)));
    assert(QOmega(0, 0, 0, 5_mpq / 64) == ring::toQOmega<ZDyadic>(ZDyadic(5, 6)));
    assert(QOmega(-7_mpq / 8, 0, 7_mpq / 8, 5_mpq / 6) ==
           ring::toQOmega<QRootTwo>(QRootTwo(5_mpq / 6, 7_mpq / 8)));
    assert(QOmega(-7_mpq / 8, 0, 7_mpq / 8, 5_mpq / 6) ==
           ring::toQOmega<QRootTwo>(QRootTwo(5_mpq / 6, 7_mpq / 8)));
    assert(QOmega(3_mpq / 512, 1_mpq / 1024, 3, 4) ==
           ring::toQOmega<DOmega>(DOmega(ZDyadic(3, 9), ZDyadic(4, 12), 3, 4)));
}

BOOST_AUTO_TEST_CASE(test_parity)
{

    assert(Z2(0) == ring::parity(1234));
    assert(Z2(1) == ring::parity(55555));
    assert(Z2(0) == ring::parity(123412341234123412341234123412341234123412341234_mpz));
    assert(Z2(1) == ring::parity(555555555555555555555555555555555555555555555555_mpz));
    assert(Z2(0) == ring::parity(ZRootTwo(1234, 1235)));
    assert(Z2(1) == ring::parity(ZRootTwo(1235, 1234)));
}

BOOST_AUTO_TEST_CASE(test_from_integer)
{
    assert(Real(2.0) == ring::fromInteger<Real>(2));
    assert(Real(2.0) == ring::fromInteger<Real>(Integer(2)));
}

template <typename T>
void test_dyadic_impl()
{
    Dyadic<T> d = Dyadic<T>(1, 5);
    Dyadic<T> d0 = Dyadic<T>(0, 9);
    Dyadic<T> dcopy = Dyadic<T>(1, 5);
    Dyadic<T> d2 = Dyadic<T>(2, 6); // 2 / 2^6 = 1 / 2^5.
    Dyadic<T> d3 = Dyadic<T>(3, 7);
    Dyadic<T> d4 = Dyadic<T>(3, 5);
    Dyadic<T> dneg = Dyadic<T>(-2, 6);

    assert(d == dcopy); // Equal with the same constructor values.
    assert(d == d2);    // Equal with different constructor values.
    assert(!(d == d3));

    Dyadic<T> d2copy = d2.copy();
    assert(d2copy == d2);
    assert(std::addressof(d2copy) != std::addressof(d2)); // Not the same object

    // assert("Dyadic(1, 5)" == d.to_string());
    // assert("Dyadic(-2, 6)" == dneg.to_string());
    //
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

    assert(Dyadic<T>(7, 7) == d + d3);

    assert(Dyadic<T>(-1, 7) == d3 - d);

    assert(Dyadic<T>(6, 13) == d2 * d3);

    assert(Dyadic<T>(-3, 7) == -d3);

    assert(1 == d.signum());
    assert(-1 == dneg.signum());
    assert(0 == d0.signum());

    assert(Dyadic<T>(1, 5) == d.abs());
    assert(Dyadic<T>(2, 6) == dneg.abs());

    assert(Dyadic<T>(1, 5) == d.adj());
    assert(Dyadic<T>(-2, 6) == dneg.adj());

    assert(Dyadic<T>(1, 5) == d.adj2());
    assert(Dyadic<T>(-2, 6) == dneg.adj2());

    Dyadic<T> dnegDenom = Dyadic<T>(3, -5);
    assert(std::make_tuple(96, 0) == dnegDenom.decompose_dyadic());
    assert(std::make_tuple(1, 5) == d.decompose_dyadic());
    assert(std::make_tuple(-1, 5) == dneg.decompose_dyadic());
    Dyadic<T> veryReducible = Dyadic<T>(128, 9);
    assert(std::make_tuple(1, 2) == veryReducible.decompose_dyadic());

    assert(2 == d.integer_of_dyadic(6));
    assert(-1 == dneg.integer_of_dyadic(5));

    assert(Dyadic<T>(6, 0) == Dyadic<T>::fromInteger(6));

    assert(Dyadic<T>(1, 5) == Dyadic<T>::fromDyadic(d));
    assert(Dyadic<T>(-2, 6) == Dyadic<T>::fromDyadic(dneg));

    assert(0.25 == ring::fromDyadic<Real>(Dyadic<T>(1, 2)));

    assert(3_mpq / 128 == ring::fromDyadic<Rational>(Dyadic<T>(3, 7)));

    assert(Dyadic<T>(1, 1) == Dyadic<T>::half());
}

BOOST_AUTO_TEST_CASE(test_dyadic)
{
    test_dyadic_impl<int>();
    test_dyadic_impl<Integer>();
}

template <typename T>
void test_roottwo_integral_impl()
{
    // rneg2 < r < r3 < r2
    RootTwo<T> r = RootTwo<T>(1, 2);
    RootTwo<T> rcopy = RootTwo<T>(1, 2);
    RootTwo<T> r2 = RootTwo<T>(4, 9);
    RootTwo<T> r3 = RootTwo<T>(-2, 5);
    RootTwo<T> rneg1 = RootTwo<T>(-3, 1);
    RootTwo<T> rneg2 = RootTwo<T>(3, -3);

    assert(r == rcopy);
    assert(!(r == r2));
    assert(r != r2);
    assert(!(r != rcopy));

    // assert("RootTwo(1, 2)" == r.to_string());
    // assert("RootTwo(4, 9)" == r2.to_string());
    //
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

    assert(RootTwo<T>(5, 2) == r + 4);

    assert(RootTwo<T>(-1, 2) == r - 2);

    assert(RootTwo<T>(12, 27) == r2 * 3);

    assert(RootTwo<T>(5, 11) == r + r2);

    assert(RootTwo<T>(-3, -7) == r - r2);

    assert(RootTwo<T>(40, 17) == r * r2);

    assert(RootTwo<T>(-1, -2) == -r);

    assert(RootTwo<T>(1, 2) == r.abs());
    assert(RootTwo<T>(-2, 5) == r3.abs());
    assert(RootTwo<T>(3, -1) == rneg1.abs());
    assert(RootTwo<T>(-3, 3) == rneg2.abs());

    assert(1 == r.signum());
    assert(1 == r3.signum());
    assert(-1 == rneg1.signum());
    assert(-1 == rneg2.signum());

    assert(RootTwo<T>(0, 1) == RootTwo<T>::roottwo());

    assert(RootTwo<T>(5, 0) == RootTwo<T>::fromInteger(5));
    assert(RootTwo<T>(-2, 0) == RootTwo<T>::fromInteger(-2));
}

BOOST_AUTO_TEST_CASE(test_roottwo_integral)
{
    test_roottwo_integral_impl<int>();
    test_roottwo_integral_impl<Integer>();
}

template <typename T>
void test_roottwo_dyadic_impl()
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

    assert(r == rcopy);
    assert(!(r == r2));
    assert(r != r2);
    assert(!(r != rcopy));

    // assert("RootTwo(Dyadic(1, 2), Dyadic(2, 3))" == r.to_string());
    // assert("RootTwo(Dyadic(1, 2), Dyadic(-2, 5))" == r2.to_string());
    //
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

    RootTwo<Dyadic<T>> sum = r + r2;
    RootTwo<Dyadic<T>> expectedSum = RootTwo<Dyadic<T>>(
        Dyadic<T>(2, 2), Dyadic<T>(6, 5));
    assert(expectedSum == sum);

    RootTwo<Dyadic<T>> difference = r - r2;
    RootTwo<Dyadic<T>> expectedDifference = RootTwo<Dyadic<T>>(
        Dyadic<T>(0, 2), Dyadic<T>(10, 5));
    assert(expectedDifference == difference);

    RootTwo<Dyadic<T>> product = r * r2;
    RootTwo<Dyadic<T>> expectedProduct = RootTwo<Dyadic<T>>(
        Dyadic<T>(8, 8), Dyadic<T>(6, 7));
    assert(expectedProduct == product);

    RootTwo<Dyadic<T>> negation = -r;
    RootTwo<Dyadic<T>> expectedNegation = RootTwo<Dyadic<T>>(
        Dyadic<T>(-1, 2), Dyadic<T>(-2, 3));
    assert(expectedNegation == negation);

    RootTwo<Dyadic<T>> expectedR0Abs = RootTwo<Dyadic<T>>(
        Dyadic<T>(0, 0), Dyadic<T>(0, 0)); // No change because r0 = 0.
    assert(expectedR0Abs == r0.abs());
    RootTwo<Dyadic<T>> expectedR2Abs = RootTwo<Dyadic<T>>(
        Dyadic<T>(1, 2), Dyadic<T>(-2, 5)); // No change because r2 is positive.
    assert(expectedR2Abs == r2.abs());
    RootTwo<Dyadic<T>> expectedR3Abs = RootTwo<Dyadic<T>>(
        Dyadic<T>(-3, 7), Dyadic<T>(2, 5)); // Get -r3 because r3 is negative.
    assert(expectedR3Abs == r3.abs());

    assert(0 == r0.signum());
    assert(1 == r.signum());
    assert(1 == r2.signum());
    assert(1 == r4.signum());
    assert(-1 == r3.signum());

    RootTwo<Dyadic<T>> expectedHalf = RootTwo<Dyadic<T>>(
        Dyadic<T>(1, 1), Dyadic<T>(0, 0));
    assert(expectedHalf == RootTwo<Dyadic<T>>::half());

    RootTwo<Dyadic<T>> expectedRootTwo = RootTwo<Dyadic<T>>(
        Dyadic<T>(0, 0), Dyadic<T>(1, 0));
    assert(expectedRootTwo == RootTwo<Dyadic<T>>::roottwo());

    RootTwo<Dyadic<T>> expectedRootHalf = RootTwo<Dyadic<T>>(
        Dyadic<T>(0, 0), Dyadic<T>(1, 1));
    assert(expectedRootHalf == RootTwo<Dyadic<T>>::roothalf());

    RootTwo<Dyadic<T>> expectedFrom5 = RootTwo<Dyadic<T>>(
        Dyadic<T>(5, 0), Dyadic<T>(0, 0));
    assert(expectedFrom5 == RootTwo<Dyadic<T>>::fromInteger(5));
    RootTwo<Dyadic<T>> expectedFromNeg2 = RootTwo<Dyadic<T>>(
        Dyadic<T>(-2, 0), Dyadic<T>(0, 0));
    assert(expectedFromNeg2 == RootTwo<Dyadic<T>>::fromInteger(-2));
}

BOOST_AUTO_TEST_CASE(test_roottwo_dyadic)
{
    test_roottwo_dyadic_impl<int>();
    test_roottwo_dyadic_impl<Integer>();
}

BOOST_AUTO_TEST_CASE(test_roottwo_rational)
{
    QRootTwo r0 = QRootTwo(0_mpq, 0_mpq);              // 0
    QRootTwo r1 = QRootTwo(2_mpq / 4, 9_mpq / 18);     // 1.207
    QRootTwo r1equal = QRootTwo(1_mpq / 2, 1_mpq / 2); // 1.207
    QRootTwo r2 = QRootTwo(-7_mpq / 13, 33_mpq / 20);  // 1.795
    QRootTwo r3 = QRootTwo(9_mpq / 17, -100_mpq / 3);  // -46.611
    QRootTwo r4 = QRootTwo(5_mpq, -7_mpq);             // -4.899

    // Make sure that we're converting the rational numbers to the standard,
    // canonical form.
    QRootTwo rNonCanonical = QRootTwo(2_mpq / 4, -3_mpq / 9);
    assert(rNonCanonical.a().get_num() == 1);
    assert(rNonCanonical.a().get_den() == 2);
    assert(rNonCanonical.b().get_num() == -1);
    assert(rNonCanonical.b().get_den() == 3);

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

    // assert("RootTwo(0, 0)" == r0.to_string());
    // assert("RootTwo(-7/13, 33/20)" == r2.to_string());
    // assert("RootTwo(9/17, -100/3)" == r3.to_string());
    // assert("RootTwo(5, -7)" == r4.to_string());
    //
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

    QRootTwo sum = r1 + r2;
    QRootTwo expectedSum = QRootTwo(-1_mpq / 26, 43_mpq / 20);
    assert(expectedSum == sum);

    QRootTwo difference = r2 - r3;
    QRootTwo expectedDifference = QRootTwo(-236_mpq / 221, 2099_mpq / 60);
    assert(expectedDifference == difference);

    QRootTwo product = r2 * r3;
    QRootTwo expectedProduct = QRootTwo(-24373_mpq / 221, 249583_mpq / 13260);
    assert(expectedProduct == product);

    QRootTwo negation = -r3;
    QRootTwo expectedNegation = QRootTwo(-9_mpq / 17, 100_mpq / 3);
    assert(expectedNegation == negation);

    QRootTwo expectedR0Abs = QRootTwo(0_mpq, 0_mpq);
    assert(QRootTwo(0_mpq, 0_mpq) == r0.abs());
    assert(QRootTwo(1_mpq / 2, 1_mpq / 2) == r1.abs());
    assert(QRootTwo(-9_mpq / 17, 100_mpq / 3) == r3.abs());

    assert(0 == r0.signum());
    assert(1 == r1.signum());
    assert(1 == r2.signum());
    assert(-1 == r3.signum());
    assert(-1 == r4.signum());

    assert(QRootTwo(1_mpq / 2, 0_mpq) == QRootTwo::half());

    assert(QRootTwo(0_mpq, 1_mpq) == QRootTwo::roottwo());

    assert(QRootTwo(0_mpq, 1_mpq / 2) == QRootTwo::roothalf());

    assert(QRootTwo(5_mpq, 0_mpq) == QRootTwo::fromInteger(5));
    assert(QRootTwo(-2_mpq, 0_mpq) == QRootTwo::fromInteger(-2));
}

BOOST_AUTO_TEST_CASE(test_zroottwo_root)
{

    assert(!ring::zroottwo_root(ZRootTwo(1, 2)).has_value());
    assert(ZRootTwo(5, 9) == ring::zroottwo_root(ZRootTwo(187, 90)));
    assert(ZRootTwo(5, 9) == ring::zroottwo_root(ZRootTwo(187, 90)));
    assert(ZRootTwo(12345678912345, 98765432198765) ==
           ring::zroottwo_root(ZRootTwo(19661636982624412467128449475_mpz, 2438652627129865854104507850_mpz)));
}

BOOST_AUTO_TEST_CASE(test_z2)
{

    assert(Z2(0) == Z2(0));
    assert(Z2(1) == Z2(1));
    assert(!(Z2(1) == Z2(0)));
    assert(!(Z2(0) == Z2(1)));
    assert(!(Z2(0) != Z2(0)));
    assert(!(Z2(1) != Z2(1)));
    assert(Z2(1) != Z2(0));
    assert(Z2(0) != Z2(1));

    // assert(Z2(0).to_string() == "Z2(0)");
    // assert(Z2(1).to_string() == "Z2(1)");
    //
    assert(Z2(0) == Z2(0) + Z2(0));
    assert(Z2(1) == Z2(0) + Z2(1));
    assert(Z2(1) == Z2(1) + Z2(0));
    assert(Z2(0) == Z2(1) + Z2(1));

    assert(Z2(0) == Z2(0) - Z2(0));
    assert(Z2(1) == Z2(0) - Z2(1));
    assert(Z2(1) == Z2(1) - Z2(0));
    assert(Z2(0) == Z2(1) - Z2(1));

    assert(Z2(0) == Z2(0) * Z2(0));
    assert(Z2(0) == Z2(0) * Z2(1));
    assert(Z2(0) == Z2(1) * Z2(0));
    assert(Z2(1) == Z2(1) * Z2(1));

    assert(-Z2(0) == Z2(0));
    assert(-Z2(1) == Z2(1));

    assert(Z2(0).abs() == Z2(0));
    assert(Z2(1).abs() == Z2(1));

    assert(Z2(0).signum() == 1);
    assert(Z2(1).signum() == 1);

    assert(Z2(0).adj() == Z2(0));
    assert(Z2(1).adj() == Z2(1));

    assert(Z2(0).adj2() == Z2(0));
    assert(Z2(1).adj2() == Z2(1));

    assert(Z2::fromInteger(0) == Z2(0));
    assert(Z2::fromInteger(1) == Z2(1));
    assert(Z2::fromInteger(123) == Z2(1));
    assert(Z2::fromInteger(222) == Z2(0));
}

template <typename T>
void test_complex_impl()
{
    Complex<T> c1 = Complex<T>(1, 2);
    Complex<T> c1copy = Complex<T>(1, 2);
    Complex<T> c2 = Complex<T>(4, 9);
    Complex<T> c3 = Complex<T>(-2, 5);
    Complex<T> c4 = Complex<T>(-3, -1);

    assert(c1 == c1copy);
    assert(!(c1 == c2));
    assert(c1 != c2);
    assert(!(c1 != c1copy));

    // if constexpr (std::is_same<T, int>::value || std::is_same<T, Integer>::value || std::is_same<T, Rational>::value)
    // {
    //     assert("Complex(1, 2)" == c1.to_string());
    //     assert("Complex(4, 9)" == c2.to_string());
    //         // }
    // else if constexpr (std::is_same<T, double>::value)
    // {
    //     assert("Complex(1.000000, 2.000000)" == c1.to_string());
    //     assert("Complex(4.000000, 9.000000)" == c2.to_string());
    //         // }
    // else if constexpr (std::is_same<T, ZDyadic>::value)
    // {
    //     assert("Complex(Dyadic(1, 0), Dyadic(2, 0))" == c1.to_string());
    //     assert("Complex(Dyadic(4, 0), Dyadic(9, 0))" == c2.to_string());
    //         // }
    // else if constexpr (std::is_same<T, DRootTwo>::value)
    // {
    //     assert("Complex(RootTwo(Dyadic(1, 0), Dyadic(0, 0)), RootTwo(Dyadic(2, 0), Dyadic(0, 0)))" == c1.to_string());
    //     assert("Complex(RootTwo(Dyadic(4, 0), Dyadic(0, 0)), RootTwo(Dyadic(9, 0), Dyadic(0, 0)))" == c2.to_string());
    //         // }
    // else if constexpr (std::is_same<T, QRootTwo>::value)
    // {
    //     assert("Complex(RootTwo(1, 0), RootTwo(2, 0))" == c1.to_string());
    //     assert("Complex(RootTwo(4, 0), RootTwo(9, 0))" == c2.to_string());
    //         // }
    // else if constexpr (std::is_same<T, Real>::value)
    // {
    // }
    // else
    // {
    //     templateError(typeid(T).name(), "Complex");
    // }

    assert(Complex<T>(5, 2) == c1 + 4);

    assert(Complex<T>(-1, 2) == c1 - 2);

    assert(Complex<T>(12, 27) == c2 * 3);

    assert(Complex<T>(5, 11) == c1 + c2);

    assert(Complex<T>(-3, -7) == c1 - c2);

    assert(Complex<T>(-53, 2) == c2 * c3);

    assert(Complex<T>(-1, -2) == -c1);

    // abs doesn't do anything for complex numbers.
    assert(Complex<T>(1, 2) == c1.abs());
    assert(Complex<T>(-2, 5) == c3.abs());

    // signum is always 1 for complex numbers.
    assert(1 == c1.signum());
    assert(1 == c2.signum());

    assert(Complex<T>(1, -2) == c1.adj());
    assert(Complex<T>(-2, -5) == c3.adj());
}

BOOST_AUTO_TEST_CASE(test_complex)
{
    test_complex_impl<int>();
    test_complex_impl<Integer>();
    test_complex_impl<Real>();
    test_complex_impl<Rational>();
    test_complex_impl<ZDyadic>();
    test_complex_impl<DRootTwo>();
    test_complex_impl<QRootTwo>();
}

template <typename T>
void test_omega_impl()
{

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

    // if constexpr (std::is_same<T, int>::value || std::is_same<T, Integer>::value || std::is_same<T, Rational>::value)
    // {
    //     assert("Omega(0, 0, 0, 0)" == o0.to_string());
    //     assert("Omega(1, -3, 3, -7)" == o1.to_string());
    //     assert("Omega(3, 7, 12, 2)" == o2.to_string());
    //     assert("Omega(5, 9, -4, -9)" == o3.to_string());
    //         // }
    // else if constexpr (std::is_same<T, double>::value)
    // {
    //     assert("Omega(0.000000, 0.000000, 0.000000, 0.000000)" == o0.to_string());
    //     assert("Omega(1.000000, -3.000000, 3.000000, -7.000000)" == o1.to_string());
    //     assert("Omega(3.000000, 7.000000, 12.000000, 2.000000)" == o2.to_string());
    //     assert("Omega(5.000000, 9.000000, -4.000000, -9.000000)" == o3.to_string());
    //         // }
    // else if constexpr (std::is_same<T, ZDyadic>::value)
    // {
    //     assert("Omega(Dyadic(0, 0), Dyadic(0, 0), Dyadic(0, 0), Dyadic(0, 0))" == o0.to_string());
    //     assert("Omega(Dyadic(1, 0), Dyadic(-3, 0), Dyadic(3, 0), Dyadic(-7, 0))" == o1.to_string());
    //     assert("Omega(Dyadic(3, 0), Dyadic(7, 0), Dyadic(12, 0), Dyadic(2, 0))" == o2.to_string());
    //     assert("Omega(Dyadic(5, 0), Dyadic(9, 0), Dyadic(-4, 0), Dyadic(-9, 0))" == o3.to_string());
    //         // }
    // else if constexpr (std::is_same<T, Real>::value)
    // {
    // }
    // else
    // {
    //     templateError(typeid(T).name(), "Omega");
    // }

    assert(Omega<T>(0, 0, 0, 0) == o0.copy());
    assert(Omega<T>(1, -3, 3, -7) == o1.copy());

    assert(Omega<T>(4, 4, 15, -5) == o1 + o2);
    assert(Omega<T>(8, 16, 8, -7) == o2 + o3);

    assert(Omega<T>(-2, -10, -9, -9) == o1 - o2);
    assert(Omega<T>(-2, -2, 16, 11) == o2 - o3);

    assert(Omega<T>(0, 0, 0, 0) == o0 * o2);
    assert(Omega<T>(-34, -22, -76, -14) == o1 * o2);
    assert(Omega<T>(63, -108, -178, -129) == o2 * o3);

    assert(Omega<T>(0, 0, 0, 0) == -o0);
    assert(Omega<T>(-1, 3, -3, 7) == -o1);
    assert(Omega<T>(-3, -7, -12, -2) == -o2);

    assert(Omega<T>(0, 0, 0, 0) == o0.abs());
    assert(Omega<T>(1, -3, 3, -7) == o1.abs());
    assert(Omega<T>(3, 7, 12, 2) == o2.abs());

    assert(1 == o0.signum());
    assert(1 == o1.signum());
    assert(1 == o2.signum());
    assert(1 == o3.signum());

    assert(Omega<T>(0, 0, 0, 0) == o0.adj());
    assert(Omega<T>(-3, 3, -1, -7) == o1.adj());
    assert(Omega<T>(-12, -7, -3, 2) == o2.adj());
    assert(Omega<T>(4, -9, -5, -9) == o3.adj());
}

BOOST_AUTO_TEST_CASE(test_omega)
{
    test_omega_impl<int>();
    test_omega_impl<Integer>();
    test_omega_impl<Real>();
    test_omega_impl<Rational>();
    test_omega_impl<ZDyadic>();
}

BOOST_AUTO_TEST_CASE(test_zroottwo_of_zomega)
{

    assert(ZRootTwo(4, -3) == ring::zroottwo_of_zomega(ZOmega(3, 0, -3, 4)));
    assert(ZRootTwo(25, 7) == ring::zroottwo_of_zomega(ZOmega(-7, 0, 7, 25)));
}

BOOST_AUTO_TEST_CASE(test_canonicalization)
{

    // Make sure Rational values are being reduced to simplest form.
    QComplex c = QComplex(555_mpq / 1665, 22_mpq / 110);
    assert(1 == c.a().get_num());
    assert(3 == c.a().get_den());
    assert(1 == c.b().get_num());
    assert(5 == c.b().get_den());
    QRootTwo r = QRootTwo(492_mpq / 1230, 66_mpq / 99);
    assert(2 == r.a().get_num());
    assert(5 == r.a().get_den());
    assert(2 == r.b().get_num());
    assert(3 == r.b().get_den());
    QOmega o = QOmega(0, 12_mpq / 24, 33_mpq / 44, 21_mpq / 49);
    assert(0 == o.a().get_num());
    assert(1 == o.a().get_den());
    assert(1 == o.b().get_num());
    assert(2 == o.b().get_den());
    assert(3 == o.c().get_num());
    assert(4 == o.c().get_den());
    assert(3 == o.d().get_num());
    assert(7 == o.d().get_den());
    QRComplex qr = QRComplex(QRootTwo(88_mpq / 99, 55_mpq / -66), QRootTwo(2, -3));
    assert(8 == qr.a().a().get_num());
    assert(9 == qr.a().a().get_den());
    assert(-5 == qr.a().b().get_num());
    assert(6 == qr.a().b().get_den());
    assert(2 == qr.b().a().get_num());
    assert(1 == qr.b().a().get_den());
    assert(-3 == qr.b().b().get_num());
    assert(1 == qr.b().b().get_den());
}
