#include "stdafx.h"
#include "regularization.h"

using namespace std;

void GetJacobian(const size_t numberSource,
	const vector<vector<vector<vector<complex<float>>>>> & a,
	const vector<vector<vector<complex<float>>>> & overline_a,
	const vector<vector<complex<float>>> & xi,
	const vector<vector<vector<complex<float>>>> & u,
	vector<vector<vector<complex<float>>>> & F_odd,
	vector<vector<vector<complex<float>>>> & F_even,
	vector<vector<complex<float>>> & F_0,
	vector<vector<complex<float>>> & F_00) noexcept
{
	size_t ii, jj;

	for (size_t count = 0; count < numberSource; ++count)
	{
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				ii = i * N + j;
				for (size_t p = 0; p < N; ++p)
				{
					for (size_t q = 0; q < N; ++q)
					{
						jj = p * N + q;
						F_odd[count][ii][jj] = a[i][j][p][q] * u[count][p][q];
						F_0[ii][jj] = a[i][j][p][q] * xi[p][q];
					}
				}
				F_0[ii][ii] += 1.0f;
			}
		}
	}

	for (size_t count = 0; count < numberSource; ++count)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					jj = p * N + q;
					F_even[count][j][jj] = overline_a[j][p][q] * u[count][p][q];
					F_00[j][jj] = overline_a[j][p][q] * xi[p][q];
				}
			}
		}
	}
}

// A_0 задачи имеет номер A[numberSource]
void GetMatrixA(const size_t numberSource,
	const vector<vector<vector<complex<float>>>> & F_odd,
	const vector<vector<vector<complex<float>>>> & F_even,
	const vector<vector<complex<float>>> & F_0,
	const vector<vector<complex<float>>> & F_00,
	vector<vector<vector<complex<float>>>> & A,
	const float alpha) noexcept
{
	vector<vector<complex<float>>> auxiliaryMatrix(N_SQUARED, vector<complex<float>>(N_SQUARED));
	MultTransposedMatrix(F_odd[0], F_odd[0], A[numberSource]);
	for (size_t count = 0; count < numberSource; ++count)
	{
		MultTransposedMatrix(F_odd[count], F_odd[count], auxiliaryMatrix);
		AddSquareMatrices(A[numberSource], auxiliaryMatrix);
	}
	for (size_t count = 0; count < numberSource; ++count)
	{
		MultTransposedMatrix(F_even[count], F_even[count], auxiliaryMatrix);
		AddSquareMatrices(A[numberSource], auxiliaryMatrix);
	}
	for (size_t count = 0; count < numberSource; ++count)
	{
		MultTransposedMatrix(F_odd[count], F_0, A[count]);
		MultTransposedMatrix(F_even[count], F_00, auxiliaryMatrix);
		AddSquareMatrices(A[count], auxiliaryMatrix);
	}

	//добавляем alpha к диагонали
	for (size_t i = 0; i < N_SQUARED; ++i)
	{
		A[numberSource][i][i] += alpha;
	}
}

void GetMatrixB(const vector<vector<complex<float>>> & F_0,
	const vector<vector<complex<float>>> & F_00,
	vector<vector<complex<float>>> & B,
	const float alpha) noexcept
{
	vector<vector<complex<float>>> auxiliaryMatrix(N_SQUARED, vector<complex<float>>(N_SQUARED));

	MultTransposedMatrix(F_0, F_0, B);
	MultTransposedMatrix(F_00, F_00, auxiliaryMatrix);
	AddSquareMatrices(B, auxiliaryMatrix);

	//добавляем alpha к диагонали
	for (size_t i = 0; i < N_SQUARED; ++i)
	{
		B[i][i] += alpha;
	}
}

void GetOperatorF(const size_t numberSource,
	const vector<vector<vector<vector<complex<float>>>>> & a,
	const vector<vector<vector<complex<float>>>> & overline_a,
	const vector<vector<complex<float>>> & xi,
	const vector<vector<vector<complex<float>>>> & u,
	const vector<vector<complex<float>>> & overline_u,
	const vector<vector<vector<complex<float>>>> & Source_R,
	const vector<vector<complex<float>>> & Source_X,
	vector<vector<complex<float>>> & F_part_odd,
	vector<vector<complex<float>>> & F_part_even) noexcept
{
	size_t ii;

	for (size_t count = 0; count < numberSource; ++count)
	{
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				ii = i * N + j;
				for (size_t p = 0; p < N; ++p)
				{
					for (size_t q = 0; q < N; ++q)
					{
						F_part_odd[count][ii] = a[i][j][p][q] * xi[p][q] * u[count][p][q];
					}
				}
				F_part_odd[count][ii] += u[count][i][j] - Source_R[count][i][j];
			}
		}
	}

	for (size_t count = 0; count < numberSource; ++count)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					F_part_even[count][j] += overline_a[j][p][q] * xi[p][q] * u[count][p][q];
				}
			}
			F_part_even[count][j] += overline_u[count][j] - Source_X[count][j];
		}
	}
}

void GetValueDerivedFunction(const size_t numberSource,
	const vector<complex<float>> & numbered_xi,
	const vector<vector<complex<float>>> & numbered_u,
	const vector<vector<vector<complex<float>>>> & F_odd,
	const vector<vector<vector<complex<float>>>> & F_even,
	const vector<vector<complex<float>>> & F_0,
	const vector<vector<complex<float>>> & F_00,
	vector<vector<complex<float>>> & F_part_odd,
	vector<vector<complex<float>>> & F_part_even) noexcept  //TODO
{
	vector<complex<float>> supportingVector(N);
	vector<complex<float>> supportingVectorSQ(N_SQUARED);

	for (size_t count = 0; count < numberSource; ++count)
	{
		MultMatrixVector(F_odd[count], numbered_xi, supportingVectorSQ);
		for (size_t i = 0; i < N_SQUARED; ++i)
		{
			F_part_odd[count][i] = supportingVectorSQ[i] - F_part_odd[count][i];
		}

		MultMatrixVector(F_0, numbered_u[count], supportingVectorSQ);
		for (size_t i = 0; i < N_SQUARED; ++i)
		{
			F_part_odd[count][i] += supportingVectorSQ[i];
		}


		MultMatrixVector(F_even[count], numbered_xi, supportingVector);
		for (size_t i = 0; i < N; ++i)
		{
			F_part_even[count][i] = supportingVector[i] - F_part_even[count][i];
		}
		MultMatrixVector(F_00, numbered_u[count], supportingVector);
		for (size_t i = 0; i < N; ++i)
		{
			F_part_even[count][i] += supportingVector[i];
		}
	}
}

// b_0 задачи находится в b[numberSource]
void Getb(const size_t numberSource,
	const vector<vector<vector<complex<float>>>> & F_odd,
	const vector<vector<vector<complex<float>>>> & F_even,
	const vector<vector<complex<float>>> & F_0,
	const vector<vector<complex<float>>> & F_00,
	const vector<vector<complex<float>>> & F_part_odd,
	const vector<vector<complex<float>>> & F_part_even,
	vector<vector<complex<float>>> & b_right) noexcept //TODO
{
	vector<complex<float>> supportingVectorSQ(N_SQUARED);

	MultTransposedMatrixVector(F_odd[0], F_part_odd[0], b_right[numberSource]);
	MultTransposedMatrixVector(F_even[0], F_part_even[0], supportingVectorSQ);
	AddVectors(b_right[numberSource], supportingVectorSQ);
	for (size_t count = 1; count < numberSource; ++count)
	{
		MultTransposedMatrixVector(F_odd[count], F_part_odd[count], supportingVectorSQ);
		AddVectors(b_right[numberSource], supportingVectorSQ);
		MultTransposedMatrixVector(F_even[count], F_part_even[count], supportingVectorSQ);
		AddVectors(b_right[numberSource], supportingVectorSQ);
	}

	for (size_t count = 0; count < numberSource; ++count)
	{
		MultTransposedMatrixVector(F_0, F_part_odd[count], b_right[count]);
		MultTransposedMatrixVector(F_00, F_part_even[count], supportingVectorSQ);
		AddVectors(b_right[count], supportingVectorSQ);
	}
}

void GetXi(const size_t numberSource,
	vector<vector<vector<complex<float>>>> & A,
	const vector<vector<complex<float>>> & inverseMatrixB,
	vector<vector<complex<float>>> & b_right,
	vector<complex<float>> & numbered_xi) noexcept
{
	vector<vector<complex<float>>> auxiliaryMatrix(N_SQUARED, vector<complex<float>>(N_SQUARED));
	vector<vector<complex<float>>> secondAuxiliaryMatrix(N_SQUARED, vector<complex<float>>(N_SQUARED));
	vector<complex<float>> supportingVectorSQ(N_SQUARED);
	vector<complex<float>> secondSupportingVectorSQ(N_SQUARED);

	//для левой части уравнения с xi все складываем в A_00 -> A[numberSource]
	for (size_t count = 0; count < numberSource; ++count)
	{
		MultMatrix(A[count], inverseMatrixB, auxiliaryMatrix);
		MultMatrixTransposed(auxiliaryMatrix, A[count], secondAuxiliaryMatrix);
		SubSquareMatrices(A[numberSource], secondAuxiliaryMatrix);
	}

	//для правой части уравнения с xi все складываем в b0 - b[numberSource]
	for (size_t count = 0; count < numberSource; ++count)
	{
		MultMatrixVector(inverseMatrixB, b_right[count], supportingVectorSQ);
		MultMatrixVector(A[count], supportingVectorSQ, secondSupportingVectorSQ);
		SubVectors(b_right[numberSource], secondSupportingVectorSQ);
	}

	// находим xi
	SolveSlauGaussa(A[numberSource], b_right[numberSource], numbered_xi);
}

void GetU(const size_t numberSource,
	const vector<vector<vector<complex<float>>>> & A,
	const vector<vector<complex<float>>> & inverseMatrixB,
	vector<vector<complex<float>>> & b_right,
	const vector<complex<float>> & numbered_xi,
	vector<vector<complex<float>>> & numbered_u) noexcept
{
	vector<complex<float>> supportingVectorSQ(N_SQUARED, complex<float>());

	for (size_t count = 0; count < numberSource; ++count)
	{
		MultTransposedMatrixVector(A[count], numbered_xi, supportingVectorSQ);
		SubVectors(b_right[count], supportingVectorSQ);
		MultMatrixVector(inverseMatrixB, b_right[count], numbered_u[count]);
	}
}
