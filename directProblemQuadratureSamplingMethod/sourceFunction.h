#pragma once
#include "basicFunctions.h"

// координаты источников, и фунция, характеризующая источник

const Point source01 = { -1.0f, 0.0f };
const Point source02 = { -1.0f, 2.5f };
const Point source03 = { -1.0f, 5.0f };
const Point source04 = { -1.0f, 7.5f };
const Point source05 = { -1.0f, 10.0f };

std::complex<float> SourceFunction(const Point source, const float x, const float y);
