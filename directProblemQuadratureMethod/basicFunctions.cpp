#include "stdafx.h"
#include "taskData.h"
#include "basicFunctions.h"


using namespace std;

float J_0(const float x)
{
	return static_cast<float>(_j0(static_cast<double>(x)));
}

float N_0(const float x)
{
	return static_cast<float>(_y0(static_cast<double>(x)));
}

complex<float> G(const Point lhs, const Point rhs)
{
	float dist = sqrt((lhs.x - rhs.x) * (lhs.x - rhs.x) + (lhs.y - rhs.y) * (lhs.y - rhs.y));
	dist *= OMEGA / C_0;
	if (dist < 0.000001f)
	{
		return static_cast<complex<float>>(0.0f);
	}
	else
	{
		return static_cast<complex<float>>(OMEGA * OMEGA * 0.25f * (N_0(dist) - I * J_0(dist)));
	}	
}

void Lasting(const string & st, clock_t & time)
{
	clock_t timeFinish = clock();
	float d = static_cast<float>(timeFinish - time) / CLOCKS_PER_SEC;
	cout << st << " " << d << endl;
	time = clock();
}
