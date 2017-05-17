#pragma once
#include "basicFunctions.h"

// координаты источников, и функция, характеризующая источник

const Posize_t source01 = { -1.0, 0.0 };
const Posize_t source02 = { -1.0, 2.5 };
const Posize_t source03 = { -1.0, 5.0 };
const Posize_t source04 = { -1.0, 7.5 };
const Posize_t source05 = { -1.0, 10.0 };

std::complex<double> SourceFunction(const Posize_t source, const double x, const double y);
