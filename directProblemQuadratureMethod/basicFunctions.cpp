#include "stdafx.h"
#include "taskData.h"
#include "basicFunctions.h"

using namespace std;

complex<float> Hankel(const float x)
{
	if (x == 0.0f)
	{
		return { 1.0f, 0.0f };
	}
	return{ static_cast<float>(_j0(x)), static_cast<float>(_y0(x)) };
}

complex<float> G(const float x_1, const float x_2, const float y_1, const float y_2)
{
	float dist = sqrt(pow(x_1 - y_1, 2) + pow(x_2 - y_2, 2));
	if(abs(dist) < 0.001f)
	{
		return{ 0.0f, 0.0f };
	}
	else
	{
		return -0.25f * I * OMEGA * OMEGA * Hankel(OMEGA * dist / C_0);
	}
}

void Lasting(const string & st, clock_t & time)
{
	clock_t timeFinish = clock();
	float d = static_cast<float>(timeFinish - time) / CLOCKS_PER_SEC;
	cout << st << " " << d << endl;
	time = clock();
}
