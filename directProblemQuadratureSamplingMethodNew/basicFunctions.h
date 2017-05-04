#pragma once

//функция Ханкеля
std::complex<float> Hankel(const float x);

// функция Грина
std::complex<float> G(const float x_1, const float x_2, const float y_1, const float y_2);

//печать времени
void Lasting(const std::string & st, const clock_t timeStart, const clock_t timeFinish);
