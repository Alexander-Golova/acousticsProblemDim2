#include "stdafx.h"
#include "Sources.h"
#include "basicFunctions.h"

using namespace std;

complex<float> Source::Function(const Point source, const float x, const float y) const
{
	float dist = sqrt(pow(x - source.x, 2) + pow(y - source.y, 2));
	return (0.0f, -0.25f) * Hankel(OMEGA * dist / C_0);
}
