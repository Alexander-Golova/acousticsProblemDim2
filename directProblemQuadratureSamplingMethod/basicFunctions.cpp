#include "stdafx.h"
#include "basicFunctions.h"

using namespace std;

complex<double> Hankel(const double x)
{
	return _j0(x) + I * _y0(x);
}

complex<double> G(const double x_1, const double x_2, const double y_1, const double y_2)
{
	double dist = sqrt(pow(x_1 - y_1, 2) + pow(x_2 - y_2, 2));
	return -0.25 * I * OMEGA * OMEGA * Hankel(OMEGA * dist / C_0);
}
