#pragma once
#include "../directProblemQuadratureSamplingMethodNew/taskData.h"
#include "../directProblemQuadratureSamplingMethod/matrix_utils.h"
#include <vector>
#include <complex>

void GetJacobian(const size_t numberSource,
	const std::vector<std::vector<std::vector<std::vector<std::complex<float>>>>> & a,
	const std::vector<std::vector<std::vector<std::complex<float>>>> & overline_a,
	const std::vector<std::vector<std::complex<float>>> & b,
	const std::vector<std::vector<std::complex<float>>> & xi,
	const std::vector<std::vector<std::vector<std::complex<float>>>> & u,
	std::vector<std::vector<std::vector<std::complex<float>>>> & F_odd,
	std::vector<std::vector<std::vector<std::complex<float>>>> & F_even);

void GetMatrixA(const size_t numberSource,
	const std::vector<std::vector<std::vector<std::complex<float>>>> & F_odd,
	const std::vector<std::vector<std::vector<std::complex<float>>>> & F_even,
	std::vector<std::vector<std::vector<std::complex<float>>>> & A, const float alpha);

void GetMatrixB(const std::vector<std::vector<std::vector<std::complex<float>>>> & F_odd,
	const std::vector<std::vector<std::vector<std::complex<float>>>> & F_even,
	std::vector<std::vector<std::complex<float>>> & B, const float alpha);

void GetOperatorF(const size_t numberSource,
	const std::vector<std::vector<std::vector<std::vector<std::complex<float>>>>> & a,
	const std::vector<std::vector<std::vector<std::complex<float>>>> & overline_a,
	const std::vector<std::vector<std::complex<float>>> & b,
	const std::vector<std::vector<std::complex<float>>> & xi,
	const std::vector<std::vector<std::vector<std::complex<float>>>> & u,
	const std::vector<std::vector<std::complex<float>>> & overline_u,
	const std::vector<std::vector<std::vector<std::complex<float>>>> & Source_R,
	const std::vector<std::vector<std::complex<float>>> & Source_X,
	std::vector<std::vector<std::complex<float>>> & F_part_odd,
	std::vector<std::vector<std::complex<float>>> & F_part_even);

void GetValueDerivedFunction(const size_t numberSource,
	const std::vector<std::complex<float>> & numbered_xi,
	const std::vector<std::vector<std::complex<float>>> & numbered_u,
	const std::vector<std::vector<std::vector<std::complex<float>>>> & F_odd,
	const std::vector<std::vector<std::vector<std::complex<float>>>> & F_even,
	std::vector<std::vector<std::complex<float>>> & F_part_odd,
	std::vector<std::vector<std::complex<float>>> & F_part_even);

void Getb(const size_t numberSource,
	const std::vector<std::vector<std::vector<std::complex<float>>>> & F_odd,
	const std::vector<std::vector<std::vector<std::complex<float>>>> & F_even,
	const std::vector<std::vector<std::complex<float>>> & F_part_odd,
	const std::vector<std::vector<std::complex<float>>> & F_part_even,
	std::vector<std::vector<std::complex<float>>> & b_right);

void GetXi(const size_t numberSource,
	std::vector<std::vector<std::vector<std::complex<float>>>> & A,
	const std::vector<std::vector<std::complex<float>>> & inverseMatrixB,
	std::vector<std::vector<std::complex<float>>> & b_right,
	std::vector<std::complex<float>> & numbered_xi);

void GetU(const size_t numberSource,
	std::vector<std::vector<std::vector<std::complex<float>>>> & A,
	const std::vector<std::vector<std::complex<float>>> & inverseMatrixB,
	std::vector<std::vector<std::complex<float>>> & b_right,
	const std::vector<std::complex<float>> & numbered_xi,
	std::vector<std::vector<std::complex<float>>> & numbered_u);