#pragma once
#include "taskData.h"
#include <vector>
#include <complex>

// источники

struct Source
{
	const size_t numberSource = 5;
	const std::vector<Point> node = {
		{ -1.0, 0.0 },{ -1.0, 2.5 },{ -1.0, 5.0 },{ -1.0, 7.5 },{ -1.0, 10.0 } };
	std::complex<double> Function(const Point source, const double x, const double y) const;
};

