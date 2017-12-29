#pragma once

#include "Sources.h"

void GetSubstantiveMatrix(const std::vector<std::vector<std::vector<std::vector<float>>>> & a,
	const std::vector<std::vector<std::vector<std::vector<float>>>> & b,
	const std::vector<std::vector<float>> & c,
	const std::vector<std::vector<float>> & xi,
	std::vector<std::vector<std::complex<float>>> & substantiveMatrix);

void GetRightPartEquation(const Source & source, const size_t count,
	std::vector<std::complex<float>> & rightPartEquation);

void InverseRenumbering(const std::vector<std::complex<float>> & numbered_u,
	std::vector<std::vector<std::complex<float>>> & u);

void GetOverlineU(const Source & source, size_t count, 
	const std::vector<std::vector<std::vector<float>>> & overline_a,
	const std::vector<std::vector<std::vector<float>>> & overline_b,
	const std::vector<std::vector<float>> & xi, const std::vector<std::vector<std::complex<float>>> & u,
	std::vector<std::complex<float>> & overline_u);
