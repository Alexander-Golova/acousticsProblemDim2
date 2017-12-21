#pragma once

const std::complex<float> I{ 0.0, 1.0 };

//функция Ханкеля
std::complex<float> Hankel(const float x);

// функция Грина
std::complex<float> G(const float x_1, const float x_2, const float y_1, const float y_2);

//печать времени
void Lasting(const std::string & st, clock_t & time);
