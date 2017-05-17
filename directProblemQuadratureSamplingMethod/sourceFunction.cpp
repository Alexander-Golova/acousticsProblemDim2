#include "stdafx.h"
#include "sourceFunction.h"

using namespace std;

complex<double> SourceFunction(const Posize_t source, const double x, const double y)
{
	double dist = sqrt(pow(x - source.x, 2) + pow(y - source.y, 2));
	return -0.25 * I * Hankel(OMEGA * dist / C_0);
}
