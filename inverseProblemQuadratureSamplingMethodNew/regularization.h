#pragma once
#include "../directProblemQuadratureSamplingMethodNew/taskData.h"
#include "../directProblemQuadratureSamplingMethod/matrix_utils.h"
#include <vector>
#include <complex>


void GetJacobian(const size_t numberSource,
	const std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> & a,
	const std::vector<std::vector<std::vector<std::complex<double>>>> & overline_a,
	const std::vector<std::vector<std::complex<double>>> & b,
	const std::vector<std::vector<std::complex<double>>> & xi,
	const std::vector<std::vector<std::vector<std::complex<double>>>> & u,
	std::vector<std::vector<std::vector<std::complex<double>>>> & F_odd,
	std::vector<std::vector<std::vector<std::complex<double>>>> & F_even);

void GetMatrixA(const size_t numberSource,
	const std::vector<std::vector<std::vector<std::complex<double>>>> & F_odd,
	const std::vector<std::vector<std::vector<std::complex<double>>>> & F_even,
	std::vector<std::vector<std::vector<std::complex<double>>>> & A, const double alpha);

void GetMatrixB(const std::vector<std::vector<std::vector<std::complex<double>>>> & F_odd,
	const std::vector<std::vector<std::vector<std::complex<double>>>> & F_even,
	std::vector<std::vector<std::complex<double>>> & B, const double alpha);

void GetOperatorF(const size_t numberSource,
	const std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> & a,
	const std::vector<std::vector<std::vector<std::complex<double>>>> & overline_a,
	const std::vector<std::vector<std::complex<double>>> & b,
	const std::vector<std::vector<std::complex<double>>> & xi,
	const std::vector<std::vector<std::vector<std::complex<double>>>> & u,
	const std::vector<std::vector<std::complex<double>>> & overline_u,
	const std::vector<std::vector<std::vector<std::complex<double>>>> & Source_R,
	const std::vector<std::vector<std::complex<double>>> & Source_X,
	std::vector<std::vector<std::complex<double>>> & F_part_odd,
	std::vector<std::vector<std::complex<double>>> & F_part_even);

void GetValueDerivedFunction(const size_t numberSource,
	const std::vector<std::complex<double>> & numbered_xi,
	const std::vector<std::vector<std::complex<double>>> & numbered_u,
	const std::vector<std::vector<std::vector<std::complex<double>>>> & F_odd,
	const std::vector<std::vector<std::vector<std::complex<double>>>> & F_even,
	std::vector<std::vector<std::complex<double>>> & F_part_odd,
	std::vector<std::vector<std::complex<double>>> & F_part_even);

void Getb(const size_t numberSource,
	const std::vector<std::vector<std::vector<std::complex<double>>>> & F_odd,
	const std::vector<std::vector<std::vector<std::complex<double>>>> & F_even,
	const std::vector<std::vector<std::complex<double>>> & F_part_odd,
	const std::vector<std::vector<std::complex<double>>> & F_part_even,
	std::vector<std::vector<std::complex<double>>> & b_right);

void GetXi(const size_t numberSource,
	std::vector<std::vector<std::vector<std::complex<double>>>> & A,
	const std::vector<std::vector<std::complex<double>>> & inverseMatrixB,
	std::vector<std::vector<std::complex<double>>> & b_right,
	std::vector<std::complex<double>> & numbered_xi);

void GetU(const size_t numberSource,
	std::vector<std::vector<std::vector<std::complex<double>>>> & A,
	const std::vector<std::vector<std::complex<double>>> & inverseMatrixB,
	std::vector<std::vector<std::complex<double>>> & b_right,
	const std::vector<std::complex<double>> & numbered_xi,
	std::vector<std::vector<std::complex<double>>> & numbered_u);
