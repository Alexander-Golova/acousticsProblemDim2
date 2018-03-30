#pragma once
#include "taskData.h"
#include <vector>
#include <complex>

//функция Бесселя нулевого порядка
float J_0(const float x);

//функция Неймана нулевого порядка
float N_0(const float x);

// функция G

std::complex<float> G(const Point lhs, const Point rhs);

//печать времени
void Lasting(const std::string & st, clock_t & time);
