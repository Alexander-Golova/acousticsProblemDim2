#include "stdafx.h"
#include "matrix_utils.h"

using namespace std;

void GetNullMatrix(vector<vector<complex<float>>> & matrix) noexcept
{
	const size_t dim1 = (size_t)matrix.size();
	const size_t dim2 = (size_t)matrix[0].size();

	for (size_t row = 0; row < dim1; ++row)
	{
		for (size_t col = 0; col < dim2; ++col)
		{
			matrix[row][col] = { 0.0f, 0.0f };
		}
	}
}

void SolveSlauGaussa(const vector<vector<complex<float>>> & matrix, const vector<complex<float>> & rhs,
	vector<complex<float>> & exactSolution) noexcept
{
	const size_t dim = (size_t)rhs.size();
	vector<vector<complex<float>>> a(dim, vector<complex<float>>(dim, complex<float>()));
	vector<complex<float>> b(dim, complex<float>());

	for (size_t i = 0; i < dim; ++i)
	{
		for (size_t j = 0; j < dim; ++j)
		{
			a[i][j] = matrix[i][j];
		}
		b[i] = rhs[i];
	}

	float maxString;
	size_t maxNomerInString;
	complex<float> c;
	complex<float> M;
	complex<float> s;

	for (size_t k = 0; k < dim; ++k)
	{
		maxString = abs(a[k][k]);
		maxNomerInString = k;

		for (size_t i = k + 1; i < dim; ++i)
		{
			if (abs(a[i][k]) > maxString)
			{
				maxString = abs(a[i][k]);
				maxNomerInString = i;
			}
		}

		swap(a[k], a[maxNomerInString]);
		swap(b[k], b[maxNomerInString]);

		for (size_t i = k + 1; i < dim; ++i)
		{
			M = a[i][k] / a[k][k];
			for (size_t j = k; j < dim; ++j)
			{
				a[i][j] -= M * a[k][j];
			}
			b[i] -= M * b[k];
		}
	}

	for (size_t i = dim - 1; i > 0; --i)
	{
		s = complex<float>();
		for (size_t j = i + 1; j < dim; ++j)
		{
			s += a[i][j] * exactSolution[j];
		}
		exactSolution[i] = (b[i] - s) / a[i][i];
	}

	s = complex<float>();
	for (size_t j = 1; j < dim; ++j)
	{
		s += a[0][j] * exactSolution[j];
	}
	exactSolution[0] = (b[0] - s) / a[0][0];
}


void InvertMatrix(vector<vector<complex<float>>> matrix, vector<vector<complex<float>>> & invertedMatrix) noexcept
{
	const size_t dim = (size_t)matrix.size();

	complex<float> temp;
	complex<float> maxElement;
	complex<float> multiplier;

	GetNullMatrix(invertedMatrix);
	for (size_t row = 0; row < dim; ++row)
	{
		invertedMatrix[row][row] = { 1.0f, 0.0f };
	}

	for (size_t k = 0; k < dim; k++)
	{
		temp = matrix[k][k];

		for (size_t j = 0; j < dim; ++j)
		{
			matrix[k][j] /= temp;
			invertedMatrix[k][j] /= temp;
		}
		for (size_t i = k + 1; i < dim; ++i)
		{
			temp = matrix[i][k];
			for (size_t j = 0; j < dim; ++j)
			{
				matrix[i][j] -= matrix[k][j] * temp;
				invertedMatrix[i][j] -= invertedMatrix[k][j] * temp;
			}
		}
	}
	for (size_t k = dim - 1; k > 0; --k)
	{
		for (size_t i = k - 1; i > 0; --i)
		{
			temp = matrix[i][k];
			for (size_t j = 0; j < dim; ++j)
			{
				matrix[i][j] -= matrix[k][j] * temp;
				invertedMatrix[i][j] -= invertedMatrix[k][j] * temp;
			}
		}
		{
			temp = matrix[0][k];
			for (size_t j = 0; j < dim; ++j)
			{
				matrix[0][j] -= matrix[k][j] * temp;
				invertedMatrix[0][j] -= invertedMatrix[k][j] * temp;
			}
		}
	}
}

void AddSquareMatrices(vector<vector<complex<float>>> & lhs, const vector<vector<complex<float>>> & rhs) noexcept
{
	const size_t dim = (size_t)lhs.size();

	for (size_t col = 0; col < dim; ++col)
	{
		for (size_t row = 0; row < dim; ++row)
		{
			lhs[col][row] += rhs[col][row];
		}
	}
}

void SubSquareMatrices(vector<vector<complex<float>>> & lhs, const vector<vector<complex<float>>> & rhs) noexcept
{
	const size_t dim = (size_t)lhs.size();

	for (size_t row = 0; row < dim; ++row)
	{
		for (size_t col = 0; col < dim; ++col)
		{
			lhs[row][col] -= rhs[row][col];
		}
	}
}

void AddVectors(vector<complex<float>> & lhs, const vector<complex<float>> & rhs) noexcept
{
	for (size_t i = 0; i < (size_t)lhs.size(); ++i)
	{
		lhs[i] += rhs[i];
	}
}

void SubVectors(vector<complex<float>> & lhs, const vector<complex<float>> & rhs) noexcept
{
	for (size_t i = 0; i < (size_t)lhs.size(); ++i)
	{
		lhs[i] -= rhs[i];
	}
}

complex<float> MultVectorVector(const vector<complex<float>> & lhs, const vector<complex<float>> & rhs) noexcept
{
	complex<float> sum = complex<float>();
	for (size_t i = 0; i < (size_t)lhs.size(); ++i)
	{
		sum += lhs[i] * rhs[i];
	}
	return sum;
}

void MultMatrixVector(const vector<vector<complex<float>>> & matrix, const vector<complex<float>> & vect,
	vector<complex<float>> & result) noexcept
{
	const size_t dim1 = (size_t)matrix.size();
	const size_t dim2 = (size_t)matrix[0].size();

	for (size_t row = 0; row < dim1; ++row)
	{
		result[row] = { 0.0f, 0.0f };
	}
	for (size_t row = 0; row < dim1; ++row)
	{
		for (size_t col = 0; col < dim2; ++col)
		{
			result[row] += matrix[row][col] * vect[col];
		}
	}
}

void MultTransposedMatrixVector(const vector<vector<complex<float>>> & matrix,
	const vector<complex<float>> & vect, vector<complex<float>> & result) noexcept
{
	const size_t dim1 = (size_t)matrix.size();
	const size_t dim2 = (size_t)matrix[0].size();

	for (size_t i = 0; i < dim2; ++i)
	{
		result[i] = { 0.0f, 0.0f };
	}
	for (size_t i = 0; i < dim2; ++i)
	{
		for (size_t j = 0; j < dim1; ++j)
		{
			result[i] += conj(matrix[j][i]) * vect[j];
		}
	}
}

void MultMatrix(const vector<vector<complex<float>>> & lhs, const vector<vector<complex<float>>> & rhs,
	vector<vector<complex<float>>> & result) noexcept
{
	const size_t dim1 = (size_t)lhs.size();
	const size_t dim2 = (size_t)lhs[0].size();
	const size_t dim3 = (size_t)rhs[0].size();

	GetNullMatrix(result);

	vector<complex<float>> thatColumn(dim3);
	vector<complex<float>> thisRow(dim2);
	complex<float> summand;

	for (size_t col = 0; col < dim3; ++col)
	{
		for (size_t inner = 0; inner < dim2; ++inner)
		{
			thatColumn[inner] = rhs[inner][col];
		}
		for (size_t row = 0; row < dim1; ++row)
		{
			thisRow = lhs[row];
			summand = { 0.0f, 0.0f };
			for (size_t inner = 0; inner < dim2; ++inner)
			{
				summand += thisRow[inner] * thatColumn[inner];
			}
			result[row][col] = summand;
		}
	}
}

void MultTransposedMatrix(const vector<vector<complex<float>>> & lhs,
	const vector<vector<complex<float>>> & rhs, vector<vector<complex<float>>> & result) noexcept
{
	const size_t dim1 = (size_t)lhs.size();
	const size_t dim2 = (size_t)lhs[0].size();
	const size_t dim3 = (size_t)rhs.size();

	vector<vector<complex<float>>> a(dim2, vector<complex<float>>(dim3, complex<float>()));

	for (size_t i = 0; i < dim2; ++i)
	{
		for (size_t j = 0; j < dim3; ++j)
		{
			a[i][j] = conj(lhs[j][i]);
		}
	}

	MultMatrix(a, rhs, result);
}

void MultMatrixTransposed(const vector<vector<complex<float>>> & lhs,
	const vector<vector<complex<float>>> & rhs, vector<vector<complex<float>>> & result) noexcept
{
	const size_t dim1 = (size_t)lhs.size();
	const size_t dim2 = (size_t)lhs[0].size();
	const size_t dim3 = (size_t)rhs.size();

	vector<vector<complex<float>>> a(dim1, vector<complex<float>>(dim3, complex<float>()));

	for (size_t i = 0; i < dim1; ++i)
	{
		for (size_t j = 0; j < dim3; ++j)
		{
			a[i][j] = conj(rhs[j][i]);
		}
	}

	MultMatrix(lhs, a, result);
}
