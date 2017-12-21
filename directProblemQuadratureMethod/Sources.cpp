#include "stdafx.h"
#include "Sources.h"
#include "basicFunctions.h"

using namespace std;

complex<float> Source::Function(const Point source, const float x, const float y) const
{
	float dist = sqrt(pow(x - source.x, 2) + pow(y - source.y, 2));
	return -0.25f * I * Hankel(OMEGA * dist / C_0);
}

void WriteSourceValues(const Source & source)
{
	ofstream fileSource("Source.txt");
	fileSource << fixed << setprecision(6);
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
		{
			for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
			{
				fileSource << source.Function(source.node[count], i * h, j * h) << " ";
			}
		}

		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			fileSource << source.Function(source.node[count], receiver, j * h) << " ";
		}
	}
	fileSource.close();
}
