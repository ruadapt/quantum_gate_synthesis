#include "../gridproblems.h"
#include <iostream>

void testLambda()
{
    std::cout << "lambda testing:" << std::endl;

    assert(ZRootTwo(1, 1) == gridprob::lambda<ZRootTwo>());
    assert(DRootTwo(1, 1) == gridprob::lambda<DRootTwo>());
    assert(QRootTwo(1, 1) == gridprob::lambda<QRootTwo>());
    std::cout << "\tlambda tests passed" << std::endl;
}

void testLambdaInv()
{
    std::cout << "lambdaInv testing:" << std::endl;

    assert(ZRootTwo(-1, 1) == gridprob::lambdaInv<ZRootTwo>());
    assert(DRootTwo(-1, 1) == gridprob::lambdaInv<DRootTwo>());
    assert(QRootTwo(-1, 1) == gridprob::lambdaInv<QRootTwo>());
    std::cout << "\tlambdaInv tests passed" << std::endl;
}

void testFloorlog()
{
    std::cout << "floorlog testing:" << std::endl;

    assert(std::make_tuple(1, 7_mpq / 5) == gridprob::floorlog<Rational>(5, 7));
    assert(std::make_tuple(0, 9) == gridprob::floorlog<Rational>(100, 9));
    Rational veryLongResult = 61174481261554107393554712930176990718174809383667600084986222987575055195730319893573610264746664394270562826004402775571158827962549405139297894619885011490074991997643815135746919626457702877321502663234254390821100188792088386046598374793226969525293787284200321130445132666864083148539066314697265625_mpq / 54232566158630092388363432844439953142084177360107193154446306624055617391387473983049776794447743968660538956361429150309481921598764136207637624928004429894199030870719225688323982375052537028496316441869000133456562964645156569800182928714876220591710828071636418845716689333746427404880400410757562368_mpz;
    assert(std::make_tuple(392, veryLongResult) == gridprob::floorlog<Rational>(6_mpq / 5, 12341234123412341234123412341234_mpz));
    std::cout << "\tfloorlog tests passed" << std::endl;
}

void testGridpointsInternal()
{
    std::cout << "gridpointsInternal testing:" << std::endl;
    assert(std::vector<ZRootTwo>{3} == gridprob::gridpointsInternal<QRootTwo>(2, 3, 3, 4));
    assert(std::vector<ZRootTwo>{3} == gridprob::gridpointsInternal<QRootTwo>(3, 3, 3, 4));
    assert(std::vector<ZRootTwo>{3} == gridprob::gridpointsInternal<QRootTwo>(2, 3, 3, 3));
    assert(std::vector<ZRootTwo>{3} == gridprob::gridpointsInternal<QRootTwo>(3, 3, 3, 3));

    {
        QRootTwo x0 = QRootTwo(1, -7);
        QRootTwo x1 = QRootTwo(2, 6);
        QRootTwo y0 = QRootTwo(100, -6);
        QRootTwo y1 = QRootTwo(200, 0);
        std::vector<ZRootTwo> points = gridprob::gridpointsInternal<QRootTwo>(x0, x1, y0, y1);
        std::sort(points.begin(), points.end());
        assert(744 == points.size());
        assert(ZRootTwo(42, -36) == points[0]);
        assert(ZRootTwo(100, -67) == points[543]);
        assert(ZRootTwo(101, -64) == points[743]);
    }

    {
        QRootTwo x0 = QRootTwo(1, -7);
        QRootTwo x1 = QRootTwo(2, 6);
        QRootTwo y0 = QRootTwo(100, -6);
        QRootTwo y1 = QRootTwo(50, 0);
        std::vector<ZRootTwo> points = gridprob::gridpointsInternal<QRootTwo>(x0, x1, y0, y1);
        std::sort(points.begin(), points.end());
        assert(std::vector<ZRootTwo>{} == points); // There are no points returned.
    }

    {
        QRootTwo x0 = QRootTwo(1, -7);
        QRootTwo x1 = QRootTwo(2, 2);
        QRootTwo y0 = QRootTwo(1, 2);
        QRootTwo y1 = QRootTwo(7, 2);
        std::vector<ZRootTwo> points = gridprob::gridpointsInternal<QRootTwo>(x0, x1, y0, y1);
        std::sort(points.begin(), points.end());
        assert(30 == points.size());
        assert(ZRootTwo(0, -6) == points[0]);
        assert(ZRootTwo(2, -2) == points[17]);
        assert(ZRootTwo(6, -1) == points[29]);
    }

    {
        QRootTwo x0 = QRootTwo(3, -6);
        QRootTwo x1 = QRootTwo(2, 2);
        QRootTwo y0 = QRootTwo(2, 4);
        QRootTwo y1 = QRootTwo(1, 5);
        std::vector<ZRootTwo> points = gridprob::gridpointsInternal<QRootTwo>(x0, x1, y0, y1);
        std::sort(points.begin(), points.end());
        assert((std::vector<ZRootTwo>{ZRootTwo(2, -4), ZRootTwo(5, -2)} == points));
    }

    {
        QRootTwo x0 = QRootTwo(4, -2);
        QRootTwo x1 = QRootTwo(2, 4);
        QRootTwo y0 = QRootTwo(2, 7);
        QRootTwo y1 = QRootTwo(15, -2);
        std::vector<ZRootTwo> points = gridprob::gridpointsInternal<QRootTwo>(x0, x1, y0, y1);
        std::sort(points.begin(), points.end());
        assert((std::vector<ZRootTwo>{ZRootTwo(5, -5), 12} == points));
    }

    {
        QRootTwo x0 = 10;
        QRootTwo x1 = 5;
        QRootTwo y0 = 9;
        QRootTwo y1 = 8;
        std::vector<ZRootTwo> points = gridprob::gridpointsInternal<QRootTwo>(x0, x1, y0, y1);
        std::sort(points.begin(), points.end());
        assert((std::vector<ZRootTwo>{} == points));
    }

    std::cout << "\tgridpointsInternal tests passed" << std::endl;
}

int main()
{
    testLambda();
    std::cout << std::endl;
    testLambdaInv();
    std::cout << std::endl;
    testFloorlog();
    std::cout << std::endl;
    testGridpointsInternal();
    std::cout << std::endl;
}