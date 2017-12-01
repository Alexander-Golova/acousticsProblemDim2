#pragma once
#include "taskData.h"
#include <vector>
#include <complex>

// источники

struct Source
{
	// количество источников
	const size_t numberSource = 5;
	// координаты источников на двухмерной плоскости
	const std::vector<Point> node = {
		{ -0.1, 0.0 }, { -0.1, 0.2 }, { -0.1, 0.4 }, { -0.1, 0.6 }, { -0.1, 0.8 }, { -0.1, 1.0 } };
	std::complex<double> Function(const Point source, const double x, const double y) const;
};

// печать значений источника в файл "Source.txt"
void WriteSourceValues(const Source & source);
