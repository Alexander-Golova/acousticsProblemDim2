#include "stdafx.h"
#include "taskData.h"
#include "basicFunctions.h"

using namespace std;

complex<double> Hankel(const double x)
{
	return{ static_cast<double>(_j0(x)), static_cast<double>(_y0(x)) };
}

complex<double> G(const double x_1, const double x_2, const double y_1, const double y_2)
{
	double dist = sqrt(pow(x_1 - y_1, 2) + pow(x_2 - y_2, 2));
	return -0.25 * I * OMEGA * OMEGA * Hankel(OMEGA * dist / C_0);
}

void Lasting(const string & st, clock_t & time)
{
	clock_t timeFinish = clock();
	double d = static_cast<double>(timeFinish - time) / CLOCKS_PER_SEC;
	cout << st << " " << d << endl;
	time = clock();
}
