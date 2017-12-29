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

void Lasting(const string & st, clock_t & time)
{
	clock_t timeFinish = clock();
	float d = static_cast<float>(timeFinish - time) / CLOCKS_PER_SEC;
	cout << st << " " << d << endl;
	time = clock();
}
