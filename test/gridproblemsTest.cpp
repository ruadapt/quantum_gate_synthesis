#define BOOST_TEST_MODULE gridProblems
#include "../gridproblems.h"
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <tuple>

namespace tt = boost::test_tools;

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
    assert(742 == points.size());
    assert(ZRootTwo(59, -48) == points.at(0));
    assert(ZRootTwo(59, -38) == points.at(543));
    assert(ZRootTwo(84, -52) == points.at(741));
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
    Operator<Integer> op = gridprob::makeOperator<Integer>(6, 7, 8, 9);
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
    Operator<Integer> op = gridprob::makeOperator<Integer>(6, 7, 8, 9);
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
    BOOST_TEST(-0.00222851183 == std::get<0>(x.value()), tt::tolerance(0.0001));
    BOOST_TEST(0.00169393714 == std::get<1>(x.value()), tt::tolerance(0.0001));
}