#pragma once

// обнуляет матрицу
void GetNullMatrix(std::vector<std::vector<std::complex<double>>> & matrix);

// решение системы методом Гаусса
void SolveSlauGaussa(const std::vector<std::vector<std::complex<double>>> & matrix,
	const std::vector<std::complex<double>> & rhs, std::vector<std::complex<double>> & exactSolution);

//нахождение обратной матрицы методом Жордана-Гаусса
void InvertMatrix(std::vector<std::vector<std::complex<double>>> matrix,
	std::vector<std::vector<std::complex<double>>> & invertedMatrix);

// сложение двух квадратных матриц - результат записывается в первую матрицу
void AddSquareMatrices(std::vector<std::vector<std::complex<double>>> & lhs,
	const std::vector<std::vector<std::complex<double>>> & rhs);

// вычитание двух квадратных матриц - результат записывается в первую матрицу
void SubSquareMatrices(std::vector<std::vector<std::complex<double>>> & lhs,
	const std::vector<std::vector<std::complex<double>>> & rhs);

// сложение двух векторов - результат записывается в первый вектор
void AddVectors(std::vector<std::complex<double>> & lhs,
	const std::vector<std::complex<double>> & rhs);

// вычитание двух векторов - результат записывается в первый вектор
void SubVectors(std::vector<std::complex<double>> & lhs,
	const std::vector<std::complex<double>> & rhs);

// умножение вектора на вектор
std::complex<double> MultVectorVector(const std::vector<std::complex<double>> & lhs,
	const std::vector<std::complex<double>> & rhs);

// умножение матрицы на вектор
void MultMatrixVector(const std::vector<std::vector<std::complex<double>>> & matrix,
	const std::vector<std::complex<double>> & vect, std::vector<std::complex<double>> & result);

// умножение транспонированной матрицы на вектор
void MultTransposedMatrixVector(const std::vector<std::vector<std::complex<double>>> & matrix,
	const std::vector<std::complex<double>> & vect, std::vector<std::complex<double>> & result);

// умножение матриц
void MultMatrix(const std::vector<std::vector<std::complex<double>>> & lhs,
				const std::vector<std::vector<std::complex<double>>> & rhs,
				std::vector<std::vector<std::complex<double>>> & result);

// умножение транспонированной матрицы на обычную матрицу
void MultTransposedMatrix(const std::vector<std::vector<std::complex<double>>> & lhs,
	const std::vector<std::vector<std::complex<double>>> & rhs,
	std::vector<std::vector<std::complex<double>>> & result);

// умножение обычной матрицы на транспонированную матрицу
void MultMatrixTransposed(const std::vector<std::vector<std::complex<double>>> & lhs,
	const std::vector<std::vector<std::complex<double>>> & rhs,
	std::vector<std::vector<std::complex<double>>> & result);
