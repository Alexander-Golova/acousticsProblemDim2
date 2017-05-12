#include "stdafx.h"
#include "sourceFunction.h"

using namespace std;

complex<float> SourceFunction(const Posize_t source, const float x, const float y)
{
	float dist = sqrt(pow(x - source.x, 2) + pow(y - source.y, 2));
	return (0.0f, -0.25f) * Hankel(OMEGA * dist / C_0);
}
