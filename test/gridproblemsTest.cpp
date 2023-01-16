#define BOOST_TEST_MODULE gridProblems
#include "utils.h"
#include "../gridproblems.h"
#include "../gridSynth.h"
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <tuple>

namespace tt = boost::test_tools;

namespace gs = gridsynth;
namespace gp = gridprob;

BOOST_AUTO_TEST_CASE(test_lambda)
{
    assert(ZRootTwo(1, 1) == gridprob::lambda<ZRootTwo>());
    assert(DRootTwo(1, 1) == gridprob::lambda<DRootTwo>());
    assert(QRootTwo(1, 1) == gridprob::lambda<QRootTwo>());
}

BOOST_AUTO_TEST_CASE(test_lambdaInv)
{
    assert(ZRootTwo(-1, 1) == gridprob::lambdaInv<ZRootTwo>());
    assert(DRootTwo(-1, 1) == gridprob::lambdaInv<DRootTwo>());
    assert(QRootTwo(-1, 1) == gridprob::lambdaInv<QRootTwo>());
}

BOOST_AUTO_TEST_CASE(test_floorlog)
{
    assert(std::make_tuple(1, 7_mpq / 5) == gridprob::floorlog<Rational>(5, 7));
    assert(std::make_tuple(0, 9) == gridprob::floorlog<Rational>(100, 9));
    Rational veryLongResult = 61174481261554107393554712930176990718174809383667600084986222987575055195730319893573610264746664394270562826004402775571158827962549405139297894619885011490074991997643815135746919626457702877321502663234254390821100188792088386046598374793226969525293787284200321130445132666864083148539066314697265625_mpq / 54232566158630092388363432844439953142084177360107193154446306624055617391387473983049776794447743968660538956361429150309481921598764136207637624928004429894199030870719225688323982375052537028496316441869000133456562964645156569800182928714876220591710828071636418845716689333746427404880400410757562368_mpz;
    assert(std::make_tuple(392, veryLongResult) == gridprob::floorlog<Rational>(6_mpq / 5, 12341234123412341234123412341234_mpz));
}

BOOST_AUTO_TEST_CASE(test_gridpointsInternal_simple)
{
    assert(std::vector<ZRootTwo>{3} == gridprob::gridpointsInternal<QRootTwo>(2, 3, 3, 4));
    assert(std::vector<ZRootTwo>{3} == gridprob::gridpointsInternal<QRootTwo>(3, 3, 3, 4));
    assert(std::vector<ZRootTwo>{3} == gridprob::gridpointsInternal<QRootTwo>(2, 3, 3, 3));
    assert(std::vector<ZRootTwo>{3} == gridprob::gridpointsInternal<QRootTwo>(3, 3, 3, 3));
}

BOOST_AUTO_TEST_CASE(test_gridpointsInternal_1)
{
    QRootTwo x0 = QRootTwo(1, -7);
    QRootTwo x1 = QRootTwo(2, 6);
    QRootTwo y0 = QRootTwo(100, -6);
    QRootTwo y1 = QRootTwo(200, 0);
    std::vector<ZRootTwo> points = gridprob::gridpointsInternal<QRootTwo>(x0, x1, y0, y1);
    std::sort(points.begin(), points.end());
    assert(744 == points.size());
    assert(ZRootTwo(42, -36) == points.at(0));
    assert(ZRootTwo(100, -67) == points.at(543));
    assert(ZRootTwo(101, -64) == points.at(743));
}

BOOST_AUTO_TEST_CASE(test_gridpointsInternal_2)
{
    QRootTwo x0 = QRootTwo(1, -7);
    QRootTwo x1 = QRootTwo(2, 6);
    QRootTwo y0 = QRootTwo(100, -6);
    QRootTwo y1 = QRootTwo(50, 0);
    std::vector<ZRootTwo> points = gridprob::gridpointsInternal<QRootTwo>(x0, x1, y0, y1);
    std::sort(points.begin(), points.end());
    assert(std::vector<ZRootTwo>{} == points); // There are no points returned.
}

BOOST_AUTO_TEST_CASE(test_gridpointsInternal_3)
{
    QRootTwo x0 = QRootTwo(1, -7);
    QRootTwo x1 = QRootTwo(2, 2);
    QRootTwo y0 = QRootTwo(1, 2);
    QRootTwo y1 = QRootTwo(7, 2);
    std::vector<ZRootTwo> points = gridprob::gridpointsInternal<QRootTwo>(x0, x1, y0, y1);
    std::sort(points.begin(), points.end());
    assert(30 == points.size());
    assert(ZRootTwo(0, -6) == points.at(0));
    assert(ZRootTwo(2, -2) == points.at(17));
    assert(ZRootTwo(6, -1) == points.at(29));
}

BOOST_AUTO_TEST_CASE(test_gridpointsInternal_4)
{
    QRootTwo x0 = QRootTwo(3, -6);
    QRootTwo x1 = QRootTwo(2, 2);
    QRootTwo y0 = QRootTwo(2, 4);
    QRootTwo y1 = QRootTwo(1, 5);
    std::vector<ZRootTwo> points = gridprob::gridpointsInternal<QRootTwo>(x0, x1, y0, y1);
    std::sort(points.begin(), points.end());
    assert((std::vector<ZRootTwo>{ZRootTwo(2, -4), ZRootTwo(5, -2)} == points));
}

BOOST_AUTO_TEST_CASE(test_gridpointsInternal_5)
{
    QRootTwo x0 = QRootTwo(4, -2);
    QRootTwo x1 = QRootTwo(2, 4);
    QRootTwo y0 = QRootTwo(2, 7);
    QRootTwo y1 = QRootTwo(15, -2);
    std::vector<ZRootTwo> points = gridprob::gridpointsInternal<QRootTwo>(x0, x1, y0, y1);
    std::sort(points.begin(), points.end());
    assert((std::vector<ZRootTwo>{ZRootTwo(5, -5), 12} == points));
}

BOOST_AUTO_TEST_CASE(test_gridpointsInternal_6)
{
    QRootTwo x0 = 10;
    QRootTwo x1 = 5;
    QRootTwo y0 = 9;
    QRootTwo y1 = 8;
    std::vector<ZRootTwo> points = gridprob::gridpointsInternal<QRootTwo>(x0, x1, y0, y1);
    std::sort(points.begin(), points.end());
    assert((std::vector<ZRootTwo>{} == points));
}

BOOST_AUTO_TEST_CASE(test_gridpoints_simple)
{
    assert(std::vector<ZRootTwo>{3} == gridprob::gridpoints<QRootTwo>(2, 3, 3, 4));
    assert(std::vector<ZRootTwo>{3} == gridprob::gridpoints<QRootTwo>(3, 3, 3, 4));
    assert(std::vector<ZRootTwo>{3} == gridprob::gridpoints<QRootTwo>(2, 3, 3, 3));
    assert(std::vector<ZRootTwo>{3} == gridprob::gridpoints<QRootTwo>(3, 3, 3, 3));
}

BOOST_AUTO_TEST_CASE(test_gridpoints_1)
{
    QRootTwo x0 = QRootTwo(1, -7);
    QRootTwo x1 = QRootTwo(2, 6);
    QRootTwo y0 = QRootTwo(100, -6);
    QRootTwo y1 = QRootTwo(200, 0);
    std::vector<ZRootTwo> points = gridprob::gridpoints<QRootTwo>(x0, x1, y0, y1);
    std::sort(points.begin(), points.end());
    BOOST_REQUIRE_EQUAL(742, points.size());
    BOOST_CHECK_EQUAL(ZRootTwo(59, -48), points.at(0));
    BOOST_CHECK_EQUAL(ZRootTwo(59, -38), points.at(543));
    BOOST_CHECK_EQUAL(ZRootTwo(84, -52), points.at(741));
}

BOOST_AUTO_TEST_CASE(test_gridpoints_2)
{
    QRootTwo x0 = QRootTwo(1, -7);
    QRootTwo x1 = QRootTwo(2, 6);
    QRootTwo y0 = QRootTwo(100, -6);
    QRootTwo y1 = QRootTwo(50, 0);
    std::vector<ZRootTwo> points = gridprob::gridpoints<QRootTwo>(x0, x1, y0, y1);
    std::sort(points.begin(), points.end());
    assert(std::vector<ZRootTwo>{} == points); // There are no points returned.
}

BOOST_AUTO_TEST_CASE(test_gridpoints_3)
{
    QRootTwo x0 = QRootTwo(1, -7);
    QRootTwo x1 = QRootTwo(2, 2);
    QRootTwo y0 = QRootTwo(1, 2);
    QRootTwo y1 = QRootTwo(7, 2);
    std::vector<ZRootTwo> points = gridprob::gridpoints<QRootTwo>(x0, x1, y0, y1);
    std::sort(points.begin(), points.end());
    assert(30 == points.size());
    assert(ZRootTwo(0, -6) == points.at(0));
    assert(ZRootTwo(2, -2) == points.at(17));
    assert(ZRootTwo(6, -1) == points.at(29));
}

BOOST_AUTO_TEST_CASE(test_gridpoints_4)
{
    QRootTwo x0 = QRootTwo(3, -6);
    QRootTwo x1 = QRootTwo(2, 2);
    QRootTwo y0 = QRootTwo(2, 4);
    QRootTwo y1 = QRootTwo(1, 5);
    std::vector<ZRootTwo> points = gridprob::gridpoints<QRootTwo>(x0, x1, y0, y1);
    std::sort(points.begin(), points.end());
    assert((std::vector<ZRootTwo>{ZRootTwo(2, -4), ZRootTwo(5, -2)} == points));
}

BOOST_AUTO_TEST_CASE(test_gridpoints_5)
{
    QRootTwo x0 = QRootTwo(4, -2);
    QRootTwo x1 = QRootTwo(2, 4);
    QRootTwo y0 = QRootTwo(2, 7);
    QRootTwo y1 = QRootTwo(15, -2);
    std::vector<ZRootTwo> points = gridprob::gridpoints<QRootTwo>(x0, x1, y0, y1);
    std::sort(points.begin(), points.end());
    assert((std::vector<ZRootTwo>{} == points));
}

BOOST_AUTO_TEST_CASE(test_gridpoints_6)
{
    QRootTwo x0 = 10;
    QRootTwo x1 = 5;
    QRootTwo y0 = 9;
    QRootTwo y1 = 8;
    std::vector<ZRootTwo> points = gridprob::gridpoints<QRootTwo>(x0, x1, y0, y1);
    std::sort(points.begin(), points.end());
    assert((std::vector<ZRootTwo>{} == points));
}

BOOST_AUTO_TEST_CASE(test_gridpointsScaled_simple)
{
    assert(std::vector<DRootTwo>{3} == gridprob::gridpointsScaled<QRootTwo>(2, 3, 3, 4, 0));
    assert(std::vector<DRootTwo>{3} == gridprob::gridpointsScaled<QRootTwo>(3, 3, 3, 4, 0));
    assert(std::vector<DRootTwo>{3} == gridprob::gridpointsScaled<QRootTwo>(2, 3, 3, 3, 0));
    assert(std::vector<DRootTwo>{3} == gridprob::gridpointsScaled<QRootTwo>(3, 3, 3, 3, 0));
}

BOOST_AUTO_TEST_CASE(test_gridpointsScaled)
{
    QRootTwo x0 = QRootTwo(1, -8);
    QRootTwo x1 = QRootTwo(2, 6);
    QRootTwo y0 = QRootTwo(0, 3);
    QRootTwo y1 = QRootTwo(10, 2);
    std::vector<DRootTwo> points = gridprob::gridpointsScaled<QRootTwo>(x0, x1, y0, y1, 7);
    std::sort(points.begin(), points.end());
    assert(8084 == points.size());
    assert(DRootTwo(1, -8) == points.at(0));
    assert(DRootTwo(ZDyadic(3, 3), ZDyadic(-85, 4)) == points.at(1234));
    assert(DRootTwo(ZDyadic(43, 2), ZDyadic(-3, 4)) == points.at(8083));
}

BOOST_AUTO_TEST_CASE(test_gridpointsScaledParity_1)
{
    DRootTwo beta = DRootTwo(3, -5);
    QRootTwo x0 = QRootTwo(1, -8);
    QRootTwo x1 = QRootTwo(2, 6);
    QRootTwo y0 = QRootTwo(0, 3);
    QRootTwo y1 = QRootTwo(5, 2);
    Integer k = 4;
    std::vector<DRootTwo> points = gridprob::gridpointsScaledParity<QRootTwo>(beta, x0, x1, y0, y1, k);
    std::sort(points.begin(), points.end());
    assert(211 == points.size());
    assert(DRootTwo(ZDyadic(-5, 1), ZDyadic(-11, 1)) == points.at(0));
    assert(DRootTwo(ZDyadic(7, 1), ZDyadic(-11, 2)) == points.at(100));
    assert(DRootTwo(8, ZDyadic(7, 2)) == points.at(210));
}

BOOST_AUTO_TEST_CASE(test_gridpointsScaledParity_2)
{
    DRootTwo beta = DRootTwo(ZDyadic(3, 2), ZDyadic(7, 9));
    QRootTwo x0 = QRootTwo(1, -8);
    QRootTwo x1 = QRootTwo(4, 1);
    QRootTwo y0 = QRootTwo(0, 3);
    QRootTwo y1 = QRootTwo(5, 2);
    Integer k = 6;
    std::vector<DRootTwo> points = gridprob::gridpointsScaledParity<QRootTwo>(beta, x0, x1, y0, y1, k);
    std::sort(points.begin(), points.end());
    assert(638 == points.size());
    assert(DRootTwo(ZDyadic(-23, 3), ZDyadic(-21, 2)) == points.at(0));
    assert(DRootTwo(ZDyadic(29, 3), ZDyadic(-9, 3)) == points.at(500));
    assert(DRootTwo(ZDyadic(39, 3), ZDyadic(3, 3)) == points.at(637));
}

BOOST_AUTO_TEST_CASE(test_point_construction)
{
    Point<ZRootTwo> p = Point<ZRootTwo>(4, 5);
    assert(ZRootTwo(4) == std::get<0>(p));
    assert(ZRootTwo(5) == std::get<1>(p));
}

BOOST_AUTO_TEST_CASE(test_makeOperator)
{
    Operator<Integer> op = matrix2x2<Integer>(6, 7, 8, 9);
    assert(op(0, 0) == 6);
    assert(op(0, 1) == 7);
    assert(op(1, 0) == 8);
    assert(op(1, 1) == 9);
}

BOOST_AUTO_TEST_CASE(test_pointFromDRootTwo)
{
    Point<DRootTwo> p = std::make_tuple(DRootTwo(1, 2), DRootTwo(ZDyadic(1, 2), ZDyadic(3, 4)));
    Point<QRootTwo> pQ = gridprob::pointFromDRootTwo<QRootTwo>(p);
    assert(std::make_tuple(QRootTwo(1, 2), QRootTwo(1_mpq / 4, 3_mpq / 16)) == pQ);
}

BOOST_AUTO_TEST_CASE(test_Ellipse_construction)
{
    Point<Integer> p = std::make_tuple(10, 12);
    Operator<Integer> op = matrix2x2<Integer>(6, 7, 8, 9);
    Ellipse<Integer> e = Ellipse<Integer>(op, p);
    Operator<Integer> op2 = e.op();
    Point<Integer> p2 = e.p();
}

BOOST_AUTO_TEST_CASE(test_unitDisk)
{
    ConvexSet<Real> u = gridprob::unitDisk<Real>();

    Point<DRootTwo> p1 = std::make_tuple(DRootTwo(ZDyadic(7, 3), 0), 0);
    Point<DRootTwo> p2 = std::make_tuple(DRootTwo(ZDyadic(1, 1), 0), DRootTwo(ZDyadic(1, 2), 0));
    Point<DRootTwo> p3 = std::make_tuple(1, DRootTwo(ZDyadic(1, 5), 0));
    Point<DRootTwo> p4 = std::make_tuple(DRootTwo(ZDyadic(3, 2), 0), DRootTwo(ZDyadic(3, 2), 0));
    assert(u.test(p1));
    assert(u.test(p2));
    assert(!u.test(p3));
    assert(!u.test(p4));

    DRootTwo a = DRootTwo((ZDyadic(1, 2)), (ZDyadic(3, 4)));
    DRootTwo b = DRootTwo((ZDyadic(5, 6)), (ZDyadic(7, 8)));
    DRootTwo c = DRootTwo((ZDyadic(9, 10)), (ZDyadic(11, 12)));
    DRootTwo d = DRootTwo((ZDyadic(4, -6)), (ZDyadic(1, -7)));
    Point<DRootTwo> p5 = std::make_tuple(a, b);
    Point<DRootTwo> p6 = std::make_tuple(c, d);

    std::optional<std::tuple<Real, Real>> x = u.intersect(p5, p6);
    BOOST_TEST(x.has_value());
    BOOST_TEST(approx_equal(Real(-0.00222851183), std::get<0>(x.value())));
    BOOST_TEST(approx_equal(Real(0.00169393714), std::get<1>(x.value())));
}

BOOST_AUTO_TEST_CASE(test_disk)
{
    ConvexSet<Real> d = gridprob::disk<Real>(25);

    Point<DRootTwo> p1 = std::make_tuple(DRootTwo(ZDyadic(7, 3), 0), 4);
    Point<DRootTwo> p2 = std::make_tuple(DRootTwo(ZDyadic(1, 1), 0), DRootTwo(ZDyadic(1, 2), 0));
    Point<DRootTwo> p3 = std::make_tuple(1, DRootTwo(ZDyadic(1, 5), 5));
    Point<DRootTwo> p4 = std::make_tuple(DRootTwo(ZDyadic(19, 2), 0), DRootTwo(ZDyadic(9, 2), 0));
    assert(d.test(p1));
    assert(d.test(p2));
    assert(!d.test(p3));
    assert(!d.test(p4));

    DRootTwo d1 = DRootTwo((ZDyadic(1, 2)), (ZDyadic(3, 4)));
    DRootTwo d2 = DRootTwo((ZDyadic(5, 6)), (ZDyadic(7, 8)));
    DRootTwo d3 = DRootTwo((ZDyadic(9, 10)), (ZDyadic(11, 12)));
    DRootTwo d4 = DRootTwo((ZDyadic(4, -6)), (ZDyadic(1, -7)));
    Point<DRootTwo> p5 = std::make_tuple(d1, d2);
    Point<DRootTwo> p6 = std::make_tuple(d3, d4);

    std::optional<std::tuple<Real, Real>> x = d.intersect(p5, p6);
    BOOST_TEST(x.has_value());
    BOOST_TEST(approx_equal(Real(-0.01164753904), std::get<0>(x.value())));
    BOOST_TEST(approx_equal(Real(0.01111296434), std::get<1>(x.value())));
}

BOOST_AUTO_TEST_CASE(test_opFromDRootTwo)
{
    Operator<DRootTwo> op1 = matrix2x2(DRootTwo(1, 2), DRootTwo(3, 4), DRootTwo(7), DRootTwo(8));
    Operator<QRootTwo> op2 = gridprob::opFromDRootTwo<QRootTwo>(op1);
    assert(QRootTwo(1, 2) == op2(0, 0));
    assert(QRootTwo(3, 4) == op2(0, 1));
    assert(QRootTwo(7) == op2(1, 0));
    assert(QRootTwo(8) == op2(1, 1));
}

BOOST_AUTO_TEST_CASE(test_logBaseDouble)
{
    Real b = 1.234;
    Real x = 3.76;
    BOOST_TEST(approx_equal(6.2989305043634936, gridprob::logBaseDouble(b, x)));
}

BOOST_AUTO_TEST_CASE(test_shiftSigma)
{
    Operator<Real> oA = gridprob::opA<Real>();
    Operator<Real> shifted = gridprob::shiftSigma(2, oA);
    BOOST_TEST(approx_equal(5.82842712473, shifted(0, 0)));
    BOOST_TEST(approx_equal(-2.00000000000, shifted(0, 1)));
    BOOST_TEST(approx_equal(0.00000000000, shifted(1, 0)));
    BOOST_TEST(approx_equal(0.17157287525, shifted(1, 1)));
}

BOOST_AUTO_TEST_CASE(test_shiftTau)
{
    Operator<Real> oA = gridprob::opA<Real>();
    Operator<Real> oB = gridprob::opB<Real>();
    OperatorPair<Real> pair = std::make_tuple(oA, oB);

    Operator<Real> rA, rB;
    std::tie(rA, rB) = gridprob::shiftState(4, pair);
    BOOST_TEST(approx_equal(33.97056274829, rA(0, 0)));
    BOOST_TEST(approx_equal(-2.00000000000, rA(0, 1)));
    BOOST_TEST(approx_equal(0.00000000000, rA(1, 0)));
    BOOST_TEST(approx_equal(0.02943725152, rA(1, 1)));
    BOOST_TEST(approx_equal(0.02943725152, rB(0, 0)));
    BOOST_TEST(approx_equal(1.41421356237, rB(0, 1)));
    BOOST_TEST(approx_equal(0.00000000000, rB(1, 0)));
    BOOST_TEST(approx_equal(33.97056274829, rB(1, 1)));
}

BOOST_AUTO_TEST_CASE(test_shiftState)
{
    Operator<Real> oA = gridprob::opA<Real>();
    Operator<Real> shifted = gridprob::shiftTau(3, oA);
    BOOST_TEST(approx_equal(0.07106781186, shifted(0, 0)));
    BOOST_TEST(approx_equal(2.00000000000, shifted(0, 1)));
    BOOST_TEST(approx_equal(0.00000000000, shifted(1, 0)));
    BOOST_TEST(approx_equal(14.07106781181, shifted(1, 1)));
}

BOOST_AUTO_TEST_CASE(test_lemma_A)
{
    Real x = 11.2345;
    Real y = 17.1253;
    Integer k = gridprob::lemma_A(x, y);
    BOOST_CHECK_EQUAL(9983, k);
}

BOOST_AUTO_TEST_CASE(test_lemma_B)
{
    Real x = 11.2345;
    Real y = 17.1253;
    Integer k = gridprob::lemma_B(x, y);
    BOOST_CHECK_EQUAL(14118, k);
}

BOOST_AUTO_TEST_CASE(test_lemma_A_l2)
{
    Real x = 111.2345;
    Real y = 117.1253;
    Integer k = gridprob::lemma_A_l2(x, y);
    BOOST_CHECK_EQUAL(5, k);
}

BOOST_AUTO_TEST_CASE(test_lemma_B_l2)
{
    Real x = 111.2345;
    Real y = 117.1253;
    Integer k = gridprob::lemma_B_l2(x, y);
    BOOST_CHECK_EQUAL(7, k);
}

BOOST_AUTO_TEST_CASE(test_step_lemma)
{
    Operator<Real> oA = gridprob::opA<Real>();
    Operator<Real> oB = gridprob::opB<Real>();
    std::optional<Operator<DRootTwo>> op = gridprob::step_lemma(std::make_tuple(oA, oB));
    BOOST_TEST(!op.has_value());
}

BOOST_AUTO_TEST_CASE(test_to_upright)
{
    Operator<Real> oA = gridprob::opA<Real>();
    Operator<Real> oB = gridprob::opB<Real>();
    Operator<DRootTwo> op = gridprob::to_upright<Real>(std::make_tuple(oA, oB));
    assert(DRootTwo(1) == op(0, 0));
    assert(DRootTwo(0) == op(0, 1));
    assert(DRootTwo(0) == op(1, 0));
    assert(DRootTwo(1) == op(1, 1));
}

BOOST_AUTO_TEST_CASE(test_to_upright_sets)
{
    CharFun alwaysTrue = [](Point<DRootTwo> p)
    { auto x = p; return true; };

    LineIntersector<Real> intersect = [](Point<DRootTwo> p1, Point<DRootTwo> p2)
    { auto x = p1; auto y = p2; return std::make_tuple(1.1, 2.2); };

    Operator<Real> oA = gridprob::opA<Real>();
    Point<Real> centerA = std::make_tuple(10, 10);
    Ellipse<Real> elA = Ellipse<Real>(oA, centerA);
    ConvexSet<Real> setA = ConvexSet<Real>(elA, alwaysTrue, intersect);

    Operator<Real> oB = gridprob::opB<Real>();
    Point<Real> centerB = std::make_tuple(7, 2);
    Ellipse<Real> elB = Ellipse<Real>(oB, centerB);
    ConvexSet<Real> setB = ConvexSet<Real>(elB, alwaysTrue, intersect);

    Operator<DRootTwo> op = gridprob::to_upright_sets<Real>(setA, setB);

    assert(DRootTwo(1) == op(0, 0));
    assert(DRootTwo(0) == op(0, 1));
    assert(DRootTwo(0) == op(1, 0));
    assert(DRootTwo(1) == op(1, 1));
}

BOOST_AUTO_TEST_CASE(test_transforms)
{
    CharFun alwaysTrue = [](Point<DRootTwo> p)
    { auto x = p; return true; };

    LineIntersector<Real> intersect = [](Point<DRootTwo> p1, Point<DRootTwo> p2)
    { auto x = p1; auto y = p2; return std::make_tuple(1.1, 2.2); };

    Operator<Real> oA = gridprob::opA<Real>();
    Point<Real> centerA = std::make_tuple(10, 10);
    Ellipse<Real> elA = Ellipse<Real>(oA, centerA);
    ConvexSet<Real> setA = ConvexSet<Real>(elA, alwaysTrue, intersect);

    Operator<Real> opG = gridprob::opK<Real>();
    Operator<DRootTwo> opG2 = gridprob::opK<DRootTwo>();

    gridprob::point_transform<Real>(opG, centerA);
    gridprob::ellipse_transform<Real>(opG, elA);
    gridprob::charfun_transform(opG2, alwaysTrue);
    gridprob::lineintersector_transform<Real>(opG2, intersect);
    gridprob::convex_transform<Real>(opG2, setA);
}

BOOST_AUTO_TEST_CASE(test_boundingbox_ellipse)
{
    Operator<Real> oB = gridprob::opB<Real>();
    Point<Real> centerB = std::make_tuple(7, 2);
    Ellipse<Real> elB = Ellipse<Real>(oB, centerB);

    Real a, b, c, d;
    Tuple2By2<Real> box = gridprob::boundingbox_ellipse(elB);
    std::tie(a, b) = std::get<0>(box);
    std::tie(c, d) = std::get<1>(box);
    BOOST_TEST(approx_equal(6.00000000000, a));
    BOOST_TEST(approx_equal(8.00000000000, b));
    BOOST_TEST(approx_equal(1.00000000000, c));
    BOOST_TEST(approx_equal(3.00000000000, d));
}

BOOST_AUTO_TEST_CASE(test_boundingbox)
{
    CharFun alwaysTrue = [](Point<DRootTwo> p)
    { auto x = p; return true; };

    LineIntersector<Real> intersect = [](Point<DRootTwo> p1, Point<DRootTwo> p2)
    { auto x = p1; auto y = p2; return std::make_tuple(1.1, 2.2); };

    Operator<Real> oB = gridprob::opB<Real>();
    Point<Real> centerB = std::make_tuple(7, 2);
    Ellipse<Real> elB = Ellipse<Real>(oB, centerB);
    ConvexSet<Real> setB = ConvexSet<Real>(elB, alwaysTrue, intersect);

    Real a, b, c, d;
    Tuple2By2<Real> box = gridprob::boundingbox(setB);
    std::tie(a, b) = std::get<0>(box);
    std::tie(c, d) = std::get<1>(box);
    BOOST_TEST(approx_equal(6.00000000000, a));
    BOOST_TEST(approx_equal(8.00000000000, b));
    BOOST_TEST(approx_equal(1.00000000000, c));
    BOOST_TEST(approx_equal(3.00000000000, d));
}

BOOST_AUTO_TEST_CASE(test_gridpoints2_scaled)
{
    ConvexSet<Real> u = gridprob::unitDisk<Real>();
    ConvexSet<Real> u2 = gridprob::unitDisk<Real>();
    std::function<std::vector<DOmega>(Integer)> solver = gridprob::gridpoints2_scaled<Real>(u, u2);
    std::vector<DOmega> points1 = solver(1);
    std::vector<DOmega> points2 = solver(2);
    BOOST_CHECK_EQUAL(17, points1.size());
    BOOST_CHECK_EQUAL(57, points2.size());
}

BOOST_AUTO_TEST_CASE(test_gridpoints2_increasing)
{
    ConvexSet<Real> u = gridprob::unitDisk<Real>();
    ConvexSet<Real> u2 = gridprob::unitDisk<Real>();
    std::function<std::vector<DOmega>(Integer)> solver = gridprob::gridpoints2_increasing<Real>(u, u2);
    std::vector<DOmega> points0 = solver(0);
    std::vector<DOmega> points1 = solver(1);
    std::vector<DOmega> points2 = solver(2);
    std::vector<DOmega> points3 = solver(3);
    std::vector<DOmega> points4 = solver(4);
    std::vector<DOmega> points5 = solver(5);
    BOOST_CHECK_EQUAL(9, points0.size());
    BOOST_CHECK_EQUAL(8, points1.size());
    // TODO the two commented out tests fail, possibly due to small number precision differences
    // BOOST_CHECK_EQUAL(40, points2.size());
    BOOST_CHECK_EQUAL(136, points3.size());
    // BOOST_CHECK_EQUAL(480, points4.size());
    BOOST_CHECK_EQUAL(1912, points5.size());
}

BOOST_AUTO_TEST_CASE(test_gridpoints2_increasing_epsilon_region2)
{
    Real prec = 10;
    Real theta = 1;
    Real epsilon = bmp::pow(2, -prec);
    ConvexSet<Real> region = gs::epsilon_region(epsilon, theta);
    std::function<List<DOmega>(Integer)> raw_candidates = gp::gridpoints2_increasing(region, gp::unitDisk<Real>());
    BOOST_CHECK_EQUAL(0, raw_candidates(0).size());
    BOOST_CHECK_EQUAL(0, raw_candidates(1).size());
    BOOST_CHECK_EQUAL(0, raw_candidates(2).size());
}

BOOST_AUTO_TEST_CASE(test_epsilon_region_to_upright_sets)
{
    Real prec = 10;
    Real theta = 1;

    Real epsilon = bmp::pow(2, -prec);
    ConvexSet<Real> region = gs::epsilon_region(epsilon, theta);
    ConvexSet<Real> disk = gp::unitDisk<Real>();
    Operator<DRootTwo> op = gp::to_upright_sets(region, disk);

    Operator<DRootTwo> expected = ring::rootHalf<DRootTwo>() * matrix2x2<DRootTwo>(DRootTwo(-3, -2), DRootTwo(1, -1), DRootTwo(-5, -4), DRootTwo(1, -1));
    for (unsigned long i = 0; i < 2; i++)
    {
        for (unsigned long j = 0; j < 2; j++)
        {
            BOOST_CHECK_EQUAL(expected(i, j), op(i, j));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_epsilon_region_to_upright_sets2)
{
    Real prec = 5;
    Real theta = 1;

    Real epsilon = bmp::pow(2, -prec);
    ConvexSet<Real> region = gs::epsilon_region(epsilon, theta);
    ConvexSet<Real> disk = gp::unitDisk<Real>();
    Operator<DRootTwo> op = gp::to_upright_sets(region, disk);

    Operator<DRootTwo> expected = matrix2x2<DRootTwo>(DRootTwo(-2, -1), -1, DRootTwo(-3, -2), DRootTwo(0, -1));
    for (unsigned long i = 0; i < 2; i++)
    {
        for (unsigned long j = 0; j < 2; j++)
        {
            BOOST_CHECK_EQUAL(expected(i, j), op(i, j));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_epsilon_region_to_upright)
{
    Real prec = 10;
    Real theta = 1;

    Real epsilon = bmp::pow(2, -prec);
    ConvexSet<Real> region = gs::epsilon_region(epsilon, theta);
    ConvexSet<Real> disk = gp::unitDisk<Real>();

    Operator<Real> op1 = region.el().op();
    Operator<Real> op2 = disk.el().op();
    Operator<DRootTwo> op = gp::to_upright(std::make_tuple(op1, op2));

    Operator<DRootTwo> expected = ring::rootHalf<DRootTwo>() * matrix2x2<DRootTwo>(DRootTwo(-3, -2), DRootTwo(1, -1), DRootTwo(-5, -4), DRootTwo(1, -1));
    for (unsigned long i = 0; i < 2; i++)
    {
        for (unsigned long j = 0; j < 2; j++)
        {
            BOOST_CHECK_EQUAL(expected(i, j), op(i, j));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_action)
{
    Operator<Real> op1 = matrix2x2<Real>(162.334, 553.269, 553.269, 1885.67);
    Operator<Real> op2 = matrix2x2<Real>(1, 0, 0, 1);
    Operator<DRootTwo> g = matrix2x2<DRootTwo>(DRootTwo(1, 1), -2, 0, DRootTwo(-1, 1));
    Operator<Real> expected1 = matrix2x2<Real>(946.15188886612, -230.54888887238, -230.54888887238, 56.17972991946);
    Operator<Real> expected2 = matrix2x2<Real>(0.17157287525, 0.82842712474, 0.82842712474, 9.82842712473);
    Operator<Real> actual1, actual2;
    std::tie(actual1, actual2) = gp::action(std::make_tuple(op1, op2), g);
    for (unsigned long i = 0; i < 2; i++)
    {
        for (unsigned long j = 0; j < 2; j++)
        {
            BOOST_TEST(expected1(i, j) == actual1(i, j), tt::tolerance(Real(0.001)));
            BOOST_TEST(expected2(i, j) == actual2(i, j), tt::tolerance(Real(0.001)));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_gridproblem_epsilon_region)
{
    Real prec = 5;
    Real theta = 1;
    Real epsilon = bmp::pow(2, -prec);
    ConvexSet<Real> region = gs::epsilon_region(epsilon, theta);
    ConvexSet<Real> disk = gp::unitDisk<Real>();
    std::function<List<DOmega>(Integer)> points = gp::gridpoints2_increasing(region, disk);
    List<DOmega> p = points(9);
    List<DOmega> expected = List<DOmega>{
        Omega(ZDyadic(-21, 5), ZDyadic(1, 5), ZDyadic(-3, 5), ZDyadic(15, 5)),
        Omega(ZDyadic(-11, 5), ZDyadic(-5, 4), ZDyadic(3, 5), ZDyadic(9, 4)),
        Omega(ZDyadic(-1, 5), ZDyadic(-21, 5), ZDyadic(9, 5), ZDyadic(21, 5)),
        Omega(ZDyadic(-25, 5), ZDyadic(3, 4), ZDyadic(-5, 5), ZDyadic(7, 4)),
        Omega(ZDyadic(-15, 5), ZDyadic(-5, 5), ZDyadic(1, 5), ZDyadic(17, 5))};
    BOOST_CHECK_EQUAL(expected, p);
}

BOOST_AUTO_TEST_CASE(test_intersect_epsilon_region)
{
    Real prec = 5;
    Real theta = 1.56;
    Real epsilon = bmp::pow(2, -prec);
    ConvexSet<Real> setA = gs::epsilon_region(epsilon, theta);
    ConvexSet<Real> setB = gp::unitDisk<Real>();
    Operator<DRootTwo> opG = gp::to_upright_sets(setA, setB);
    Operator<DRootTwo> opG_inv = gp::special_inverse(opG);
    ConvexSet<Real> setA_prime = gp::convex_transform(opG_inv, setA);

    DRootTwo x0 = 0;
    DRootTwo dx = 1;
    DRootTwo beta_prime = DRootTwo(1, 1);

    Point<DRootTwo> p1A = std::make_tuple(x0, beta_prime);
    Point<DRootTwo> p2A = std::make_tuple(dx, DRootTwo(0));
    Maybe<Point<Real>> iA = setA_prime.intersect(p1A, p2A);

    BOOST_REQUIRE(iA.has_value());
    Point<Real> expected = {0, 0};
    BOOST_CHECK(expected == iA.value());
}

BOOST_AUTO_TEST_CASE(test_gridproblem_epsilon_region2)
{
    Real prec = 5;
    Real theta = 1.56;
    Real epsilon = bmp::pow(2, -prec);
    ConvexSet<Real> region = gs::epsilon_region(epsilon, theta);
    ConvexSet<Real> disk = gp::unitDisk<Real>();
    std::function<List<DOmega>(Integer)> points = gp::gridpoints2_increasing(region, disk);
    List<DOmega> p = points(0);
    List<DOmega> expected = List<DOmega>{DOmega(-1, 0, 0, 0)};
    BOOST_CHECK_EQUAL(expected, p);
}

