#pragma once

#include "../directProblemQuadratureSamplingMethodNew/taskData.h"
#include <vector>
#include <complex>

void ProjectionXi(std::vector<std::vector<std::complex<float>>> & xi);

void PrintXi(std::vector<std::vector<std::complex<float>>> & xi, size_t iteration);
