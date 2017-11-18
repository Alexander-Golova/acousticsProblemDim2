#include "stdafx.h"
#include "Sources.h"
#include "basicFunctions.h"

using namespace std;

complex<double> Source::Function(const Point source, const double x, const double y) const
{
	double dist = sqrt(pow(x - source.x, 2) + pow(y - source.y, 2));
	return -0.25 * I * Hankel(OMEGA * dist / C_0);
}

void WriteSourceValues(const Source & source)
{
	ofstream fileSource("Source.txt");
	fileSource << fixed << setprecision(6);
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		for (size_t i = 0; i <= SPLITTING; ++i)
		{
			for (size_t j = 0; j <= SPLITTING; ++j)
			{
				fileSource << source.Function(source.node[count], i * h, j * h) << " ";
			}
		}

		for (size_t i = 0; i <= SPLITTING; ++i)
		{
			for (size_t j = 0; j <= SPLITTING; ++j)
			{
				fileSource << source.Function(source.node[count], RECEIVER + i * stepReceiver, j * h) << " ";
			}
		}
	}
	fileSource.close();
}
