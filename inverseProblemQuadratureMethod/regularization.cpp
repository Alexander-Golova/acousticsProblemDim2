#include "stdafx.h"
#include "regularization.h"

using namespace std;

void GetJacobian(const size_t numberSource, const vector<vector<vector<vector<complex<double>>>>> & a,
	const vector<vector<vector<complex<double>>>> & overline_a, const vector<vector<complex<double>>> & b,
	const vector<vector<complex<double>>> & xi, const vector<vector<vector<complex<double>>>> & u,
	vector<vector<vector<complex<double>>>> & F_odd, vector<vector<vector<complex<double>>>> & F_even)
{
	complex<double> sumOfTheCoefficients;
	size_t ii, jj;

	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			ii = i * (NUMBER_PARTITION_POINT + 1) + j;
			sumOfTheCoefficients = { 0.0, 0.0 };
			for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
			{
				for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
				{
					jj = p * (NUMBER_PARTITION_POINT + 1) + q;
					if ((i != p) || (q != j))
					{
						sumOfTheCoefficients += a[i][j][p][q];
						for (size_t count = 1; count <= numberSource; ++count)
						{
							F_odd[count][ii][jj] = a[i][j][p][q] * u[count - 1][p][q];
						}
						F_odd[0][ii][jj] = a[i][j][p][q] * xi[p][q];
					}
				}
			}
			for (size_t count = 1; count <= numberSource; ++count)
			{
				F_odd[count][ii][ii] = u[count - 1][i][j] * (b[i][j] - sumOfTheCoefficients);
			}
			F_odd[0][ii][ii] = 1.0;
			F_odd[0][ii][ii] += xi[i][j] * (b[i][j] - sumOfTheCoefficients);
		}
	}

	for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
	{
		for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
		{
			for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
			{
				jj = p * (NUMBER_PARTITION_POINT + 1) + q;
				for (size_t count = 1; count <= numberSource; ++count)
				{
					F_even[count][j][jj] = overline_a[j][p][q] * u[count - 1][p][q];
				}
				F_even[0][j][jj] = overline_a[j][p][q] * xi[p][q];
			}
		}
	}
}

void GetMatrixA(const size_t numberSource,
	const vector<vector<vector<complex<double>>>> & F_odd, const vector<vector<vector<complex<double>>>> & F_even,
	vector<vector<vector<complex<double>>>> & A, const double alpha)
{
	const size_t N_squared = (NUMBER_PARTITION_POINT + 1) * (NUMBER_PARTITION_POINT + 1);
	vector<vector<complex<double>>> auxiliaryMatrix(N_squared, vector<complex<double>>(N_squared, complex<double>()));
	MultTransposedMatrix(F_odd[1], F_odd[1], A[0]);
	for (size_t count = 2; count <= numberSource; ++count)
	{
		MultTransposedMatrix(F_odd[count], F_odd[count], auxiliaryMatrix);
		AddSquareMatrices(A[0], auxiliaryMatrix);
	}
	for (size_t count = 1; count <= numberSource; ++count)
	{
		MultTransposedMatrix(F_even[count], F_even[count], auxiliaryMatrix);
		AddSquareMatrices(A[0], auxiliaryMatrix);
	}
	for (size_t count = 1; count <= numberSource; ++count)
	{
		MultTransposedMatrix(F_odd[count], F_odd[0], A[count]);
		MultTransposedMatrix(F_even[count], F_even[0], auxiliaryMatrix);
		AddSquareMatrices(A[count], auxiliaryMatrix);
	}

	//добавляем alpha к диагонали
	for (size_t i = 0; i < N_squared; ++i)
	{
		A[0][i][i] += alpha;
	}
}

void GetMatrixB(const vector<vector<vector<complex<double>>>> & F_odd,
	const vector<vector<vector<complex<double>>>> & F_even,
	vector<vector<complex<double>>> & B, const double alpha)
{
	const size_t N_squared = (NUMBER_PARTITION_POINT + 1) * (NUMBER_PARTITION_POINT + 1);
	vector<vector<complex<double>>> auxiliaryMatrix(N_squared, vector<complex<double>>(N_squared, complex<double>()));

	MultTransposedMatrix(F_odd[0], F_odd[0], B);
	MultTransposedMatrix(F_even[0], F_even[0], auxiliaryMatrix);
	AddSquareMatrices(B, auxiliaryMatrix);

	//добавляем alpha к диагонали
	for (size_t i = 0; i < N_squared; ++i)
	{
		B[i][i] += alpha;
	}
}

void GetOperatorF(const size_t numberSource, const vector<vector<vector<vector<complex<double>>>>> & a,
	const vector<vector<vector<complex<double>>>> & overline_a, const vector<vector<complex<double>>> & b,
	const vector<vector<complex<double>>> & xi, const vector<vector<vector<complex<double>>>> & u,
	const vector<vector<complex<double>>> & overline_u, const vector<vector<vector<complex<double>>>> & Source_R,
	const vector<vector<complex<double>>> & Source_X, vector<vector<complex<double>>> & F_part_odd,
	vector<vector<complex<double>>> & F_part_even)
{
	size_t ii;
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			ii = i * (NUMBER_PARTITION_POINT + 1) + j;
			for (size_t count = 0; count < numberSource; ++count)
			{
				F_part_odd[count][ii] = u[count][i][j];
			}
			for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
			{
				for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
				{
					if ((i != p) || (q != j))
					{
						for (size_t count = 0; count < numberSource; ++count)
						{
							F_part_odd[count][ii] += a[i][j][p][q] * (xi[p][q] * u[count][p][q] - xi[i][j] * u[count][i][j]);
						}
					}
				}
			}
			for (size_t count = 0; count < numberSource; ++count)
			{
				F_part_odd[count][ii] += b[i][j] * xi[i][j] * u[count][i][j];
				F_part_odd[count][ii] -= Source_R[count][i][j];
			}
		}
	}
	for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
	{
		for (size_t count = 0; count < numberSource; ++count)
		{
			F_part_even[count][j] = overline_u[count][j] - Source_X[count][j];
		}
		for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
		{
			for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
			{
				for (size_t count = 0; count < numberSource; ++count)
				{
					F_part_even[count][j] += overline_a[j][p][q] * xi[p][q] * u[count][p][q];
				}
			}
		}
	}
}

void GetValueDerivedFunction(const size_t numberSource, const vector<complex<double>> & numbered_xi,
	const vector<vector<complex<double>>> & numbered_u, const vector<vector<vector<complex<double>>>> & F_odd,
	const vector<vector<vector<complex<double>>>> & F_even, vector<vector<complex<double>>> & F_part_odd,
	vector<vector<complex<double>>> & F_part_even)
{
	const size_t N_squared = (NUMBER_PARTITION_POINT + 1) * (NUMBER_PARTITION_POINT + 1);
	vector<complex<double>> supportingVector_square(N_squared, complex<double>());
	vector<complex<double>> supportingVector(N_squared, complex<double>());

	for (size_t count = 1; count <= numberSource; ++count)
	{
		MultMatrixVector(F_odd[count], numbered_xi, supportingVector_square);
		for (size_t i = 0; i < N_squared; ++i)
		{
			F_part_odd[count - 1][i] = supportingVector_square[i] - F_part_odd[count - 1][i];
		}
		MultMatrixVector(F_odd[0], numbered_u[count - 1], supportingVector_square);
		for (size_t i = 0; i < N_squared; ++i)
		{
			F_part_odd[count - 1][i] += supportingVector_square[i];
		}
		MultMatrixVector(F_even[count], numbered_xi, supportingVector);

		for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
		{
			F_part_even[count - 1][i] = supportingVector_square[i] - F_part_even[count - 1][i];
		}
		MultMatrixVector(F_even[0], numbered_u[count - 1], supportingVector);
		for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
		{
			F_part_even[count - 1][i] += supportingVector_square[i];
		}
	}
}

void Getb(const size_t numberSource, const vector<vector<vector<complex<double>>>> & F_odd,
	const vector<vector<vector<complex<double>>>> & F_even, const vector<vector<complex<double>>> & F_part_odd,
	const vector<vector<complex<double>>> & F_part_even, vector<vector<complex<double>>> & b_right)
{
	const size_t N_squared = (NUMBER_PARTITION_POINT + 1) * (NUMBER_PARTITION_POINT + 1);
	vector<complex<double>> supportingVector_square(N_squared, complex<double>());

	MultTransposedMatrixVector(F_odd[1], F_part_odd[0], b_right[0]);
	MultTransposedMatrixVector(F_even[1], F_part_even[0], supportingVector_square);
	AddVectors(b_right[0], supportingVector_square);
	for (size_t count = 2; count <= numberSource; ++count)
	{
		MultTransposedMatrixVector(F_odd[count], F_part_odd[count - 1], supportingVector_square);
		AddVectors(b_right[0], supportingVector_square);
		MultTransposedMatrixVector(F_even[count], F_part_even[count - 1], supportingVector_square);
		AddVectors(b_right[0], supportingVector_square);
	}

	for (size_t count = 1; count <= numberSource; ++count)
	{
		MultTransposedMatrixVector(F_odd[0], F_part_odd[count - 1], b_right[count]);
		MultTransposedMatrixVector(F_even[0], F_part_even[count - 1], supportingVector_square);
		AddVectors(b_right[count], supportingVector_square);
	}
}

void GetXi(const size_t numberSource, vector<vector<vector<complex<double>>>> & A,
	const vector<vector<complex<double>>> & inverseMatrixB, vector<vector<complex<double>>> & b_right,
	vector<complex<double>> & numbered_xi)
{
	const size_t N_squared = (NUMBER_PARTITION_POINT + 1) * (NUMBER_PARTITION_POINT + 1);
	vector<vector<complex<double>>> auxiliaryMatrix(N_squared, vector<complex<double>>(N_squared, complex<double>()));
	vector<vector<complex<double>>> secondAuxiliaryMatrix(N_squared, vector<complex<double>>(N_squared, complex<double>()));
	vector<complex<double>> supportingVector_square(N_squared, complex<double>());
	vector<complex<double>> secondSupportingVector_square(N_squared, complex<double>());

	//для левой части уравнения с xi все складываем в A_00
	for (size_t count = 1; count <= numberSource; ++count)
	{
		MultMatrix(A[count], inverseMatrixB, auxiliaryMatrix);
		MultMatrixTransposed(auxiliaryMatrix, A[count], secondAuxiliaryMatrix);
		SubSquareMatrices(A[0], secondAuxiliaryMatrix);
	}

	//для правой части уравнения с xi все складываем в b0
	for (size_t count = 1; count <= numberSource; ++count)
	{
		MultMatrixVector(inverseMatrixB, b_right[count], supportingVector_square);
		MultMatrixVector(A[count], supportingVector_square, secondSupportingVector_square);
		SubVectors(b_right[0], secondSupportingVector_square);
	}

	// находим xi
	SolveSlauGaussa(A[0], b_right[0], numbered_xi);
}

void GetU(const size_t numberSource, vector<vector<vector<complex<double>>>> & A,
	const vector<vector<complex<double>>> & inverseMatrixB, vector<vector<complex<double>>> & b_right,
	const vector<complex<double>> & numbered_xi, vector<vector<complex<double>>> & numbered_u)
{
	const size_t N_squared = (NUMBER_PARTITION_POINT + 1) * (NUMBER_PARTITION_POINT + 1);
	vector<complex<double>> supportingVector_square(N_squared, complex<double>());

	for (size_t count = 0; count < numberSource; ++count)
	{
		MultTransposedMatrixVector(A[count + 1], numbered_xi, supportingVector_square);
		SubVectors(b_right[count + 1], supportingVector_square);
		MultMatrixVector(inverseMatrixB, b_right[count + 1], numbered_u[count]);
	}
}
