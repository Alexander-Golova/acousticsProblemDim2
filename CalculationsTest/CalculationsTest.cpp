#include "stdafx.h"
#include "..\directProblemQuadratureMethod\basicFunctions.h"

using namespace std;

BOOST_AUTO_TEST_SUITE(Hankel_function)

// Вычисляем функцию Ханкеля
BOOST_AUTO_TEST_CASE(computes_the_Hankel_function)
{
	complex<float> value = { 1.0f, 0.0f };
	BOOST_CHECK(value == Hankel(0.0f));
}
BOOST_AUTO_TEST_SUITE_END(Hankel_function)
