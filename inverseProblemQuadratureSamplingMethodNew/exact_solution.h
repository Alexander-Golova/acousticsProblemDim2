#pragma once

#include "../directProblemQuadratureSamplingMethodNew/taskData.h"
#include <vector>
#include <complex>

void ProjectionXi(std::vector<std::vector<std::complex<float>>> & xi);

void PrintXi(std::vector<std::vector<std::complex<float>>> & xi, size_t iteration);

void Renumbering(const std::vector<std::vector<std::complex<float>>> & xi, std::vector<std::complex<float>> & numbered_xi);

void InverseRenumbering(const std::vector<std::complex<float>> & numbered_xi, std::vector<std::vector<std::complex<float>>> & xi);
