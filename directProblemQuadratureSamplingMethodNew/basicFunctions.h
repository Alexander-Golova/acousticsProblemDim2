#pragma once

const std::complex<double> I(0.0, 1.0);

//������� �������
std::complex<double> Hankel(const double x);

// ������� �����
std::complex<double> G(const double x_1, const double x_2, const double y_1, const double y_2);

//������ �������
void Lasting(const std::string & st, const clock_t timeStart, const clock_t timeFinish);
