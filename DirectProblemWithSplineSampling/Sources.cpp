#include "stdafx.h"
#include "Sources.h"
#include "basicFunctions.h"

using namespace std;

complex<double> Source::Function(const Point source, const double x, const double y) const
{
	double dist = sqrt(pow(x - source.x, 2) + pow(y - source.y, 2));
	return -0.25 * I * Hankel(OMEGA * dist / C_0);
}
