#pragma once

//������� �������
std::complex<float> Hankel(const float x);

// ������� �����
std::complex<float> G(const float x_1, const float x_2, const float y_1, const float y_2);

//������ �������
void Lasting(const std::string & st, const clock_t timeStart, const clock_t timeFinish);
