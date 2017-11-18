#pragma once

void GetExactSolution(std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1> & xi);

void WriteSolutionFile(std::array<std::array<std::complex<double>, SPLITTING + 1>, SPLITTING + 1> & xi);
