#pragma once

#include "../directProblemQuadratureMethod/taskData.h"
#include <vector>
#include <complex>

void ProjectionXi(std::vector<std::vector<std::complex<float>>> & xi);

void PrintXi(std::vector<std::vector<std::complex<float>>> & xi, size_t iteration);

void RenumberingXi(const std::vector<std::vector<std::complex<float>>> & xi, std::vector<std::complex<float>> & numbered_xi);

void RenumberingU(const std::vector<std::vector<std::complex<float>>> & u, std::vector<std::complex<float>> & numbered_u);

void InverseRenumberingXi(const std::vector<std::complex<float>> & numbered_xi, std::vector<std::vector<std::complex<float>>> & xi);

void InverseRenumberingU(const std::vector<std::complex<float>> & numbered_xi, std::vector<std::vector<std::complex<float>>> & xi);
