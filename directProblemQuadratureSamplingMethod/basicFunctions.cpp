#include "stdafx.h"
#include "basicFunctions.h"

using namespace std;

complex<float> Hankel(const float x)
{
	return (static_cast<float>(_j0(x)), static_cast<float>(_y0(x)));
}

complex<float> G(const float x_1, const float x_2, const float y_1, const float y_2)
{
	float dist = sqrt(pow(x_1 - y_1, 2) + pow(x_2 - y_2, 2));
	return (0.0f, -0.25f * OMEGA * OMEGA) * Hankel(OMEGA * dist / C_0);
}
