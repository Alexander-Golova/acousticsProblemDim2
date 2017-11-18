#pragma once

#include <complex>
#include <array>
#include "basicArrays.h"
#include "Sources.h"

void GetSubstantiveMatrix(BasicArrays & basicArrays, std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1> & xi,
	std::array<std::array<std::complex<double>, N_SQUARED>, N_SQUARED> & substantiveMatrix);

void GetRightPartEquation(const Source & source, size_t count,
	std::array<std::complex<double>, N_SQUARED> & rightPartEquation);

void InverseRenumbering(std::array<std::complex<double>, N_SQUARED> & numbered_u, std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1> & u);

void GetOverlineU(const Source & source, size_t count, BasicArrays & basicArrays,
	std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1> & xi,
	std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1> & u,
	std::array<std::array<std::complex<double>, N_SQUARED>, N_SQUARED> & overline_u);
