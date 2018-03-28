#pragma once
#include "taskData.h"
#include <vector>
#include <complex>

void GetBasicArrays(std::vector<std::vector<std::vector<std::vector<float>>>> & a,
	std::vector<std::vector<std::vector<std::vector<float>>>> & b,
	std::vector<std::vector<float>> & c,
	std::vector<std::vector<std::vector<float>>> & overline_a,
	std::vector<std::vector<std::vector<float>>> & overline_b) noexcept;
