#pragma once
#include "basicFunctions.h"

// координаты источников, и фунция, характеризующая источник

const Posize_t source01 = { -1.0f, 0.0f };
const Posize_t source02 = { -1.0f, 2.5f };
const Posize_t source03 = { -1.0f, 5.0f };
const Posize_t source04 = { -1.0f, 7.5f };
const Posize_t source05 = { -1.0f, 10.0f };

std::complex<float> SourceFunction(const Posize_t source, const float x, const float y);
