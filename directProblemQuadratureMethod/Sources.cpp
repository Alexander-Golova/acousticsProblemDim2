#include "stdafx.h"
#include "Sources.h"
#include "basicFunctions.h"

using namespace std;

complex<float> Source::Function(const Point source, const float x, const float y) const
{
	const float dist = OMEGA * sqrt(static_cast<float>((x - source.x) * (x - source.x) + (y - source.y) * (y - source.y))) / C_0;

	return 0.25f * static_cast<complex<float>>(N_0(dist) - I * J_0(dist));
}

void WriteSourceValues(const Source & source)
{
	ofstream fileSource("Source.txt");
	fileSource << fixed << setprecision(6);
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				fileSource << source.Function(source.node[count], i * step, j * step) << " ";
			}
		}

		for (size_t j = 0; j < N; ++j)
		{
			fileSource << source.Function(source.node[count], receiver, j * step) << " ";
		}
	}
	fileSource.close();
}
