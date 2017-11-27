#pragma once

#include "Sources.h"

void GetSubstantiveMatrix(const std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> & a,
	const std::vector<std::vector<std::complex<double>>> & b,
	const std::vector<std::vector<double>> & xi,
	std::vector<std::vector<std::complex<double>>> & substantiveMatrix);

void GetRightPartEquation(const Source & source, size_t count,
	std::vector<std::complex<double>> & rightPartEquation);

void InverseRenumbering(const std::vector<std::complex<double>> & numbered_u,
	std::vector<std::vector<std::complex<double>>> & u);

void GetOverlineU(const Source & source, size_t count, 
	const std::vector<std::vector<std::vector<std::complex<double>>>> & overline_a,
	const std::vector<std::vector<double>> & xi, const std::vector<std::vector<std::complex<double>>> & u,
	std::vector<std::complex<double>> & overline_u);
