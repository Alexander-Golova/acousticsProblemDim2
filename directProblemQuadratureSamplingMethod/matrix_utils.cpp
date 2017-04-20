#include "stdafx.h"
#include "matrix_utils.h"

using namespace std;


void GetNullMatrix(vector<vector<complex<float>>> & matrix)
{
	const size_t dim1 = (size_t)matrix.size();
	const size_t dim2 = (size_t)matrix[0].size();

	for (size_t row = 0; row < dim1; ++row)
	{
		for (size_t col = 0; col < dim2; ++col)
		{
			matrix[row][col] = (0.0f, 0.0f);
		}
	}
}

void SolveSlauGaussa(const vector<vector<complex<float>>> & matrix, const vector<complex<float>> & rhs,
						vector<complex<float>> & exactSolution)
{
	const size_t dim = (size_t)rhs.size();
	vector<vector<complex<float>>> a(dim, vector<complex<float>>(dim + 1, complex<float>()));

	for (size_t row = 0; row < dim; ++row)
	{
		for (size_t col = 0; col < dim; ++col)
		{
			a[row][col] = matrix[row][col];
		}
		a[row][dim] = rhs[row];
	}

	vector<size_t> where(dim);
	size_t sel;

	for (size_t col = 0, row = 0; col < dim && row <= dim; ++col)
	{
		sel = row;
		for (size_t i = row; i < dim; ++i)
		{
			if (abs(a[i][col]) > abs(a[sel][col]))
			{
				sel = i;
			}
		}

		swap(a[sel], a[row]);
		/*for (size_t i = col; i <= dim; ++i)
		{
			swap(a[sel][i], a[row][i]);
		}*/
		where[col] = row;

		for (size_t i = 0; i < dim; ++i)
		{
			if (i != row)
			{
				complex<float> c = a[i][col] / a[row][col];
				for (size_t j = col; j <= dim; ++j)
				{
					a[i][j] -= a[row][j] * c;
				}
			}
		}
		++row;
	}

	exactSolution.assign(dim, complex<float>());
	for (size_t i = 0; i < dim; ++i)
	{
		exactSolution[i] = a[where[i]][dim] / a[where[i]][i];
	}
}

void InvertMatrix(vector<vector<complex<float>>> matrix, vector<vector<complex<float>>> & invertedMatrix)
{
	const size_t dim = (size_t)matrix.size();

	complex<float> temp;
	complex<float> maxElement;
	complex<float> multiplier;

	size_t maxNumber;

	GetNullMatrix(invertedMatrix);
	for (size_t row = 0; row < dim; ++row)
	{
		invertedMatrix[row][row] = (1.0f, 0.0f);
	}

	for (size_t row = 0; row < dim; ++row)
	{
		maxElement = matrix[row][row];
		maxNumber = row;
		for (size_t col = row; col < dim; ++col)
		{
			if (abs(matrix[col][row]) > abs(maxElement))
			{
				maxElement = matrix[col][row];
				maxNumber = col;
			}
		}
		for (size_t col = 0; col < dim; ++col)
		{
			temp = matrix[row][col];
			matrix[row][col] = matrix[maxNumber][col];
			matrix[maxNumber][col] = temp;
		}
		for (size_t col = 0; col < dim; ++col)
		{
			temp = invertedMatrix[row][col];
			invertedMatrix[row][col] = invertedMatrix[maxNumber][col];
			invertedMatrix[maxNumber][col] = temp;
		}

		for (size_t col = 0; col < dim; ++col)
		{
			if (row != col)
			{
				multiplier = matrix[col][row] / matrix[row][row];
				for (size_t i = 0; i < dim; ++i)
				{
					matrix[col][i] -= matrix[row][i] * multiplier;
					invertedMatrix[col][i] -= invertedMatrix[row][i] * multiplier;
				}
			}
		}
	}
	for (size_t row = 0; row < dim; ++row)
	{
		invertedMatrix[row][row] /= matrix[row][row];
	}
}

void AddSquareMatrices(vector<vector<complex<float>>> & lhs, const vector<vector<complex<float>>> & rhs)
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

void SubSquareMatrices(vector<vector<complex<float>>> & lhs, const vector<vector<complex<float>>> & rhs)
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

void AddVectors(vector<complex<float>> & lhs, const vector<complex<float>> & rhs)
{
	for (size_t i = 0; i < (size_t)lhs.size(); ++i)
	{
		lhs[i] += rhs[i];
	}
}

void SubVectors(vector<complex<float>> & lhs, const vector<complex<float>> & rhs)
{
	for (size_t i = 0; i < (size_t)lhs.size(); ++i)
	{
		lhs[i] -= rhs[i];
	}
}

complex<float> MultVectorVector(const vector<complex<float>> & lhs, const vector<complex<float>> & rhs)
{
	complex<float> sum = complex<float>();
	for (size_t i = 0; i < (size_t)lhs.size(); ++i)
	{
		sum += lhs[i] * rhs[i];
	}
	return sum;
}

void MultMatrixVector(const vector<vector<complex<float>>> & matrix, const vector<complex<float>> & vect,
	vector<complex<float>> & result)
{
	const size_t dim1 = (size_t)matrix.size();
	const size_t dim2 = (size_t)matrix[0].size();

	for (size_t row = 0; row < dim1; ++row)
	{
		result[row] = (0.0f, 0.0f);
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
	const vector<complex<float>> & vect, vector<complex<float>> & result)
{
	const size_t dim1 = (size_t)matrix.size();
	const size_t dim2 = (size_t)matrix[0].size();

	for (size_t i = 0; i < dim2; ++i)
	{
		result[i] = (0.0f, 0.0f);
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
	vector<vector<complex<float>>> & result)
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
			summand = (0.0f, 0.0f);
			for (size_t inner = 0; inner < dim2; ++inner)
			{
				summand += thisRow[inner] * thatColumn[inner];
			}
			result[row][col] = summand;
		}
	}
}

void MultTransposedMatrix(const vector<vector<complex<float>>> & lhs,
	const vector<vector<complex<float>>> & rhs, vector<vector<complex<float>>> & result)
{
	const size_t dim1 = (size_t)lhs.size();
	const size_t dim2 = (size_t)lhs[0].size();
	const size_t dim3 = (size_t)rhs.size();

	GetNullMatrix(result);

	for (size_t row = 0; row < dim2; ++row)
	{
		for (size_t col = 0; col < dim3; ++col)
		{
			for (size_t inner = 0; inner < dim1; ++inner)
			{
				result[row][col] += conj(lhs[inner][row]) * rhs[inner][col];
			}
		}
	}
}

void MultMatrixTransposed(const vector<vector<complex<float>>> & lhs,
	const vector<vector<complex<float>>> & rhs, vector<vector<complex<float>>> & result)
{
	const size_t dim1 = (size_t)lhs.size();
	const size_t dim2 = (size_t)lhs[0].size();
	const size_t dim3 = (size_t)rhs.size();

	GetNullMatrix(result);

	for (size_t row = 0; row < dim1; ++row)
	{
		for (size_t col = 0; col < dim3; ++col)
		{
			for (size_t inner = 0; inner < dim2; ++inner)
			{
				result[row][col] += lhs[row][inner] * conj(rhs[inner][col]);
			}
		}
	}
}
