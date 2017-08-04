#pragma once
#include "taskData.h"
#include <vector>
#include <complex>

// ���������

struct Source
{
	// ���������� ����������
	const size_t numberSource = 5;
	// ���������� ���������� �� ���������� ���������
	const std::vector<Point> node = {
		{ -0.1, 0.0 }, { -0.1, 2.0 }, { -0.1, 4.0 }, { -0.1, 6.0 }, { -0.1, 8.0 }, { -0.1, 10.0 } };
	std::complex<double> Function(const Point source, const double x, const double y) const;
};

// ������ �������� ��������� � ���� "Source.txt"
void WriteSourceValues(const Source & source);
