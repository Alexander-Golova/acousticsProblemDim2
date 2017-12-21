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
		{ -0.1f, 0.0f }, { -0.1f, 0.2f }, { -0.1f, 0.4f }, { -0.1f, 0.6f }, { -0.1f, 0.8f }, { -0.1f, 1.0f } };
	std::complex<float> Function(const Point source, const float x, const float y) const;
};

// печать значений источника в файл "Source.txt"
void WriteSourceValues(const Source & source);
