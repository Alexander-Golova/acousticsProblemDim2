#pragma once

#include "../directProblemQuadratureSamplingMethodNew/taskData.h"
#include <vector>
#include <complex>

void ProjectionXi(std::vector<std::vector<std::complex<double>>> & xi);

void PrintXi(std::vector<std::vector<std::complex<double>>> & xi, size_t iteration);

void Renumbering(const std::vector<std::vector<std::complex<double>>> & xi, std::vector<std::complex<double>> & numbered_xi);

void InverseRenumbering(const std::vector<std::complex<double>> & numbered_xi, std::vector<std::vector<std::complex<double>>> & xi);
