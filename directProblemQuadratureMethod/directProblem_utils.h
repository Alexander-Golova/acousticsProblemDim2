#pragma once

#include "Sources.h"

void GetSubstantiveMatrix(const std::vector<std::vector<std::vector<std::vector<std::complex<float>>>>> & a,
	const std::vector<std::vector<float>> & xi,
	std::vector<std::vector<std::complex<float>>> & substantiveMatrix) noexcept;

void GetRightPartEquation(const Source & source, const size_t count,
	std::vector<std::complex<float>> & rightPartEquation) noexcept;

void InverseRenumbering(const std::vector<std::complex<float>> & numbered_u,
	std::vector<std::vector<std::complex<float>>> & u) noexcept;

void GetOverlineU(const Source & source, size_t count, 
	const std::vector<std::vector<std::vector<std::complex<float>>>> & overline_a,
	const std::vector<std::vector<float>> & xi, const std::vector<std::vector<std::complex<float>>> & u,
	std::vector<std::complex<float>> & overline_u) noexcept;
