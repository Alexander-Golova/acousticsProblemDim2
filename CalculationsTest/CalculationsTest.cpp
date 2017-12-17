#include "stdafx.h"
#include "..\directProblemQuadratureMethod\basicFunctions.h"

using namespace std;

BOOST_AUTO_TEST_SUITE(Hankel_function)

// Вычисляем функцию Ханкеля
BOOST_AUTO_TEST_CASE(computes_the_Hankel_function)
{
	complex<double> value = { 1.0, 0.0 };
	cout << Hankel(0.0);
	BOOST_CHECK(value == Hankel(0.0));
}

BOOST_AUTO_TEST_SUITE_END()
