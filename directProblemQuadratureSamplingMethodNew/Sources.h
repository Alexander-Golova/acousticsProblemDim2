#pragma once
#include "taskData.h"
#include <vector>
#include <complex>

// источники

struct Source
{
	const size_t numberSource = 5;
	const std::vector<Point> node = {
		{ -1.0f, 0.0f }, { -1.0f, 2.5f }, { -1.0f, 5.0f }, { -1.0f, 7.5f }, { -1.0f, 10.0f } };
	std::complex<float> Function(const Point source, const float x, const float y) const;
};
