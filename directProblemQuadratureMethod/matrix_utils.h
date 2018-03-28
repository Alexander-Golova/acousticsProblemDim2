#pragma once

// обнуляет матрицу
void GetNullMatrix(std::vector<std::vector<std::complex<float>>> & matrix) noexcept;

// решение системы методом Гаусса
void SolveSlauGaussa(const std::vector<std::vector<std::complex<float>>> & matrix,
	const std::vector<std::complex<float>> & rhs, std::vector<std::complex<float>> & exactSolution) noexcept;

//нахождение обратной матрицы методом Жордана-Гаусса
void InvertMatrix(std::vector<std::vector<std::complex<float>>> matrix,
	std::vector<std::vector<std::complex<float>>> & invertedMatrix) noexcept;

// сложение двух квадратных матриц - результат записывается в первую матрицу
void AddSquareMatrices(std::vector<std::vector<std::complex<float>>> & lhs,
	const std::vector<std::vector<std::complex<float>>> & rhs) noexcept;

// вычитание двух квадратных матриц - результат записывается в первую матрицу
void SubSquareMatrices(std::vector<std::vector<std::complex<float>>> & lhs,
	const std::vector<std::vector<std::complex<float>>> & rhs) noexcept;

// сложение двух векторов - результат записывается в первый вектор
void AddVectors(std::vector<std::complex<float>> & lhs,
	const std::vector<std::complex<float>> & rhs) noexcept;

// вычитание двух векторов - результат записывается в первый вектор
void SubVectors(std::vector<std::complex<float>> & lhs,
	const std::vector<std::complex<float>> & rhs) noexcept;

// умножение вектора на вектор
std::complex<float> MultVectorVector(const std::vector<std::complex<float>> & lhs,
	const std::vector<std::complex<float>> & rhs) noexcept;

// умножение матрицы на вектор
void MultMatrixVector(const std::vector<std::vector<std::complex<float>>> & matrix,
	const std::vector<std::complex<float>> & vect, std::vector<std::complex<float>> & result) noexcept;

// умножение транспонированной матрицы на вектор
void MultTransposedMatrixVector(const std::vector<std::vector<std::complex<float>>> & matrix,
	const std::vector<std::complex<float>> & vect, std::vector<std::complex<float>> & result) noexcept;

// умножение матриц
void MultMatrix(const std::vector<std::vector<std::complex<float>>> & lhs,
				const std::vector<std::vector<std::complex<float>>> & rhs,
				std::vector<std::vector<std::complex<float>>> & result) noexcept;

// умножение транспонированной матрицы на обычную матрицу
void MultTransposedMatrix(const std::vector<std::vector<std::complex<float>>> & lhs,
	const std::vector<std::vector<std::complex<float>>> & rhs,
	std::vector<std::vector<std::complex<float>>> & result) noexcept;

// умножение обычной матрицы на транспонированную матрицу
void MultMatrixTransposed(const std::vector<std::vector<std::complex<float>>> & lhs,
	const std::vector<std::vector<std::complex<float>>> & rhs,
	std::vector<std::vector<std::complex<float>>> & result) noexcept;
