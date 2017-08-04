#pragma once

#include "taskData.h"
#include <vector>
#include <complex>
#include <array>

struct BasicArrays
{
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> aa_1;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> ab_1;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> ac_1;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> bb_1;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> bc_1;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> cc_1;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> overline_aa_1;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> overline_ab_1;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> overline_ac_1;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> overline_bb_1;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> overline_bc_1;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> overline_cc_1;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> aa_2;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> ab_2;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> ac_2;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> bb_2;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> bc_2;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> cc_2;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> overline_aa_2;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> overline_ab_2;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> overline_ac_2;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> overline_bb_2;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> overline_bc_2;
	std::array<std::array<std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1>, SPLITTING + 1> overline_cc_2;
};

void GetBasicArrays(BasicArrays basicArrays);

void WriteBasicArraysFile(BasicArrays basicArrays);

