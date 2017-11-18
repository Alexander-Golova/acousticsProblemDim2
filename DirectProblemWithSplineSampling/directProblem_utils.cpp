#include "stdafx.h"
#include "taskData.h"
#include "directProblem_utils.h"
#include "basicArrays.h"
#include "Sources.h"

using namespace std;

void GetSubstantiveMatrix(BasicArrays & basicArrays, array<array<complex<double>, SPLITTING + 1>, SPLITTING + 1> & xi,
	array<array<complex<double>, N_SQUARED>, N_SQUARED> & substantiveMatrix)
{
	size_t ii, jj;
	size_t auxInd;
	const size_t N = SPLITTING;

	// считаем для вершин квадрата
	for (size_t i = 0; i <= SPLITTING; ++i)
	{
		for (size_t j = 0; j <= SPLITTING; ++j)
		{
			ii = i * (SPLITTING + 1) + j;
			
			// вершина O(0; 0)
			//1 треугольник
			substantiveMatrix[ii][0] = basicArrays.aa_1[i][j][0][0] * xi[0][0];
			substantiveMatrix[ii][0] += basicArrays.ab_1[i][j][0][0] * xi[0][1];
			substantiveMatrix[ii][0] += basicArrays.ac_1[i][j][0][0] * xi[1][0];

			//2 треугольник
			// вершина A(N; 0)
			auxInd = SPLITTING * (SPLITTING + 1);
			//1 треугольник
			substantiveMatrix[ii][auxInd] = basicArrays.aa_1[i][j][N][0] * xi[N][0];
			substantiveMatrix[ii][auxInd] += basicArrays.ab_1[i][j][N][0] * xi[N][1];
			substantiveMatrix[ii][auxInd] += basicArrays.ac_1[i][j][N - 1][0] * xi[N - 1][0];
			substantiveMatrix[ii][auxInd] += basicArrays.bc_1[i][j][N - 1][0] * xi[N - 1][1];
			substantiveMatrix[ii][auxInd] += basicArrays.cc_1[i][j][N - 1][0] * xi[N][0];
			//2 треугольник
			substantiveMatrix[ii][auxInd] += basicArrays.ac_2[i][j][N - 1][0] * xi[N][1];
			substantiveMatrix[ii][auxInd] += basicArrays.bc_2[i][j][N - 1][0] * xi[N - 1][1];
			substantiveMatrix[ii][auxInd] += basicArrays.cc_2[i][j][N - 1][0] * xi[N][0];

			// вершина B(0; N)
			auxInd = SPLITTING;
			//1 треугольник
			substantiveMatrix[ii][auxInd] = basicArrays.aa_1[i][j][0][N] * xi[0][N];
			substantiveMatrix[ii][auxInd] += basicArrays.ab_1[i][j][0][N - 1] * xi[0][N - 1];
			substantiveMatrix[ii][auxInd] += basicArrays.ac_1[i][j][0][N] * xi[1][N];
			substantiveMatrix[ii][auxInd] += basicArrays.bb_1[i][j][0][N - 1] * xi[0][N];
			substantiveMatrix[ii][auxInd] += basicArrays.bc_1[i][j][0][N - 1] * xi[1][N - 1];
			//2 треугольник
			substantiveMatrix[ii][auxInd] += basicArrays.ab_2[i][j][0][N - 1] * xi[1][N];
			substantiveMatrix[ii][auxInd] += basicArrays.bb_2[i][j][0][N - 1] * xi[0][N];
			substantiveMatrix[ii][auxInd] += basicArrays.bc_2[i][j][0][N - 1] * xi[1][N - 1];

			// вершина C(N; N)
			auxInd = SPLITTING * (SPLITTING + 1) + SPLITTING;
			//1 треугольник
			substantiveMatrix[ii][auxInd] = basicArrays.aa_1[i][j][N][N] * xi[N][N];
			substantiveMatrix[ii][auxInd] += basicArrays.ac_1[i][j][N - 1][N] * xi[N - 1][N];
			substantiveMatrix[ii][auxInd] += basicArrays.bb_1[i][j][N][N - 1] * xi[N][N];
			substantiveMatrix[ii][auxInd] += basicArrays.cc_1[i][j][N - 1][N] * xi[N][N];
			//2 треугольник
			substantiveMatrix[ii][auxInd] += basicArrays.aa_2[i][j][N - 1][N - 1] * xi[N][N];
			substantiveMatrix[ii][auxInd] += basicArrays.ab_2[i][j][N - 1][N] * xi[N][N];
			substantiveMatrix[ii][auxInd] += basicArrays.ac_2[i][j][N - 1][N - 1] * xi[N][N - 1];
			substantiveMatrix[ii][auxInd] += basicArrays.bb_2[i][j][N][N - 1] * xi[N][N];
			substantiveMatrix[ii][auxInd] += basicArrays.cc_2[i][j][N - 1][N] * xi[N][N];
		}
	}

	// считаем для сторон квадрата
	for (size_t i = 0; i <= SPLITTING; ++i)
	{
		for (size_t j = 0; j <= SPLITTING; ++j)
		{
			ii = i * (SPLITTING + 1) + j;

			// сторона OA(p; 0)
			for (size_t p = 1; p < SPLITTING; ++p)
			{
				auxInd = p * (SPLITTING + 1);
				//1 треугольник
				substantiveMatrix[ii][auxInd] = basicArrays.aa_1[i][j][p][0] * xi[p][0];
				substantiveMatrix[ii][auxInd] += basicArrays.ab_1[i][j][p][0] * xi[p][1];
				substantiveMatrix[ii][auxInd] += basicArrays.ac_1[i][j][p - 1][0] * xi[p - 1][0];
				substantiveMatrix[ii][auxInd] += basicArrays.ac_1[i][j][p][0] * xi[p + 1][0];
				substantiveMatrix[ii][auxInd] += basicArrays.bc_1[i][j][p - 1][0] * xi[p - 1][1];
				substantiveMatrix[ii][auxInd] += basicArrays.cc_1[i][j][p - 1][0] * xi[p][0];
				//2 треугольник
				substantiveMatrix[ii][auxInd] += basicArrays.ac_2[i][j][p - 1][0] * xi[p][1];
				substantiveMatrix[ii][auxInd] += basicArrays.bc_2[i][j][p - 1][0] * xi[p - 1][1];
				substantiveMatrix[ii][auxInd] += basicArrays.cc_2[i][j][p - 1][0] * xi[p][0];
			}

			// сторона OB(0; s)
			for (size_t s = 1; s < SPLITTING; ++s)
			{
				auxInd = s;
				//1 треугольник
				substantiveMatrix[ii][auxInd] = basicArrays.aa_1[i][j][0][s] * xi[0][s];
				substantiveMatrix[ii][auxInd] += basicArrays.ab_1[i][j][0][s - 1] * xi[0][s - 1];
				substantiveMatrix[ii][auxInd] += basicArrays.ab_1[i][j][0][s] * xi[0][s + 1];
				substantiveMatrix[ii][auxInd] += basicArrays.ac_1[i][j][0][s] * xi[1][s];
				substantiveMatrix[ii][auxInd] += basicArrays.bb_1[i][j][0][s - 1] * xi[0][s];
				substantiveMatrix[ii][auxInd] += basicArrays.bc_1[i][j][0][s - 1] * xi[1][s - 1];
				//2 треугольник
				substantiveMatrix[ii][auxInd] += basicArrays.ab_2[i][j][0][s - 1] * xi[1][s];
				substantiveMatrix[ii][auxInd] += basicArrays.bb_2[i][j][0][s - 1] * xi[0][s];
				substantiveMatrix[ii][auxInd] += basicArrays.bc_2[i][j][0][s - 1] * xi[1][s - 1];
			}
			//
			// сторона BC(p; N)
			//
			for (size_t p = 1; p < SPLITTING; ++p)
			{
				auxInd = p * (SPLITTING + 1) + SPLITTING;
				//1 треугольник
				substantiveMatrix[ii][auxInd] = basicArrays.aa_1[i][j][p][N] * xi[p][N];
				substantiveMatrix[ii][auxInd] += basicArrays.ac_1[i][j][p - 1][N] * xi[p - 1][N];
				substantiveMatrix[ii][auxInd] += basicArrays.ac_1[i][j][p][N] * xi[p + 1][N];
				substantiveMatrix[ii][auxInd] += basicArrays.bb_1[i][j][p][N - 1] * xi[p][N];
				substantiveMatrix[ii][auxInd] += basicArrays.bc_1[i][j][p][N - 1] * xi[p + 1][N - 1];
				substantiveMatrix[ii][auxInd] += basicArrays.cc_1[i][j][p - 1][N] * xi[p][N];
				//2 треугольник
				substantiveMatrix[ii][auxInd] += basicArrays.aa_2[i][j][p - 1][N - 1] * xi[p][N];
				substantiveMatrix[ii][auxInd] += basicArrays.ab_2[i][j][p][N - 1] * xi[p + 1][N];
				substantiveMatrix[ii][auxInd] += basicArrays.ab_2[i][j][p - 1][N - 1] * xi[p - 1][N];
				substantiveMatrix[ii][auxInd] += basicArrays.ac_2[i][j][p - 1][N - 1] * xi[p][N - 1];
				substantiveMatrix[ii][auxInd] += basicArrays.bb_2[i][j][p][N - 1] * xi[p][N];
				substantiveMatrix[ii][auxInd] += basicArrays.bc_2[i][j][p][N - 1] * xi[p + 1][N - 1];
				substantiveMatrix[ii][auxInd] += basicArrays.cc_2[i][j][p - 1][N] * xi[p][N];
			}

			// сторона AC(N; s)
			for (size_t s = 1; s < SPLITTING; ++s)
			{
				auxInd = SPLITTING * (SPLITTING + 1) + s;
				//1 треугольник
				substantiveMatrix[ii][auxInd] = basicArrays.aa_1[i][j][N][s] * xi[N][s];
				substantiveMatrix[ii][auxInd] += basicArrays.ab_1[i][j][N][s - 1] * xi[N][s - 1];
				substantiveMatrix[ii][auxInd] += basicArrays.ab_1[i][j][N][s] * xi[N][s + 1];
				substantiveMatrix[ii][auxInd] += basicArrays.ac_1[i][j][N - 1][s] * xi[N - 1][s];
				substantiveMatrix[ii][auxInd] += basicArrays.bb_1[i][j][N][s - 1] * xi[N][s];
				substantiveMatrix[ii][auxInd] += basicArrays.bc_1[i][j][N - 1][s] * xi[N - 1][s + 1];
				substantiveMatrix[ii][auxInd] += basicArrays.cc_1[i][j][N - 1][s] * xi[N][s];
				//2 треугольник
				substantiveMatrix[ii][auxInd] += basicArrays.aa_2[i][j][N - 1][s - 1] * xi[N][s];
				substantiveMatrix[ii][auxInd] += basicArrays.ab_2[i][j][N - 1][s - 1] * xi[N - 1][s];
				substantiveMatrix[ii][auxInd] += basicArrays.ac_2[i][j][N - 1][s] * xi[N][s + 1];
				substantiveMatrix[ii][auxInd] += basicArrays.ac_2[i][j][N - 1][s - 1] * xi[N][s - 1];
				substantiveMatrix[ii][auxInd] += basicArrays.bb_2[i][j][N][s - 1] * xi[N][s];
				substantiveMatrix[ii][auxInd] += basicArrays.bc_2[i][j][N - 1][s] * xi[N - 1][s + 1];
				substantiveMatrix[ii][auxInd] += basicArrays.cc_2[i][j][N - 1][s] * xi[N][s];
			}
		}
	}

	// считаем для внутренних точек квадрата
	for (size_t i = 1; i < SPLITTING; ++i)
	{
		for (size_t j = 1; j < SPLITTING; ++j)
		{
			ii = i * (SPLITTING + 1) + j;
			for (size_t p = 1; p < SPLITTING; ++p)
			{
				for (size_t s = 1; s < SPLITTING; ++s)
				{
					jj = p * (SPLITTING + 1) + s;
					//1 треугольник
					substantiveMatrix[ii][jj] = basicArrays.aa_1[i][j][p][s] * xi[p][s];
					substantiveMatrix[ii][jj] += basicArrays.ab_1[i][j][p][s - 1] * xi[p][s - 1];
					substantiveMatrix[ii][jj] += basicArrays.ab_1[i][j][p][s] * xi[p][s + 1];
					substantiveMatrix[ii][jj] += basicArrays.ac_1[i][j][p - 1][s] * xi[p - 1][s];
					substantiveMatrix[ii][jj] += basicArrays.ac_1[i][j][p][s] * xi[p + 1][s];
					substantiveMatrix[ii][jj] += basicArrays.bb_1[i][j][p][s - 1] * xi[p][s];
					substantiveMatrix[ii][jj] += basicArrays.bc_1[i][j][p - 1][s] * xi[p - 1][s + 1];
					substantiveMatrix[ii][jj] += basicArrays.bc_1[i][j][p][s - 1] * xi[p + 1][s - 1];
					substantiveMatrix[ii][jj] += basicArrays.cc_1[i][j][p - 1][s] * xi[p][s];
					//2 треугольник
					substantiveMatrix[ii][jj] += basicArrays.aa_2[i][j][p - 1][s - 1] * xi[p][s];
					substantiveMatrix[ii][jj] += basicArrays.ab_2[i][j][p][s - 1] * xi[p + 1][s];
					substantiveMatrix[ii][jj] += basicArrays.ab_2[i][j][p - 1][s - 1] * xi[p - 1][s];
					substantiveMatrix[ii][jj] += basicArrays.ac_2[i][j][p - 1][s] * xi[p][s + 1];
					substantiveMatrix[ii][jj] += basicArrays.ac_2[i][j][p - 1][s - 1] * xi[p][s - 1];
					substantiveMatrix[ii][jj] += basicArrays.bb_2[i][j][p][s - 1] * xi[p][s];
					substantiveMatrix[ii][jj] += basicArrays.bc_2[i][j][p - 1][s] * xi[p - 1][s + 1];
					substantiveMatrix[ii][jj] += basicArrays.bc_2[i][j][p][s - 1] * xi[p + 1][s - 1];
					substantiveMatrix[ii][jj] += basicArrays.cc_2[i][j][p - 1][s] * xi[p][s];
				}
			}
		}
	}
	for (size_t i = 0; i < N_SQUARED; ++i)
	{
		substantiveMatrix[i][i] += 1.0;
	}
}

void GetRightPartEquation(const Source & source, size_t count, array<complex<double>, N_SQUARED> & rightPartEquation)
{
	size_t ii;
	for (size_t i = 0; i <= SPLITTING; ++i)
	{
		for (size_t j = 0; j <= SPLITTING; ++j)
		{
			ii = i * (SPLITTING + 1) + j;
			rightPartEquation[ii] = source.Function(source.node[count], i * h, j * h);
		}
	}
}

void InverseRenumbering(array<complex<double>, N_SQUARED> & numbered_u, array<array<complex<double>, SPLITTING + 1>, SPLITTING + 1> & u)
{
	size_t coordinate_x;
	size_t coordinate_y;
	for (size_t i = 0; i < N_SQUARED; ++i)
	{
		coordinate_x = i / (SPLITTING + 1);
		coordinate_y = i % (SPLITTING + 1);
		u[coordinate_x][coordinate_y] = numbered_u[i];
	}
}

void GetOverlineU(const Source & source, size_t count, BasicArrays & basicArrays,
	array<array<complex<double>, SPLITTING + 1>, SPLITTING + 1> & xi,
	array<array<complex<double>, SPLITTING + 1>, SPLITTING + 1> & u,
	array<array<complex<double>, N_SQUARED>, N_SQUARED> & overline_u)
{
	for (size_t i = 0; i <= SPLITTING; ++i)
	{
		for (size_t j = 0; j <= SPLITTING; ++j)
		{
			overline_u[i][j] = source.Function(source.node[count], RECEIVER + i * stepReceiver, j * h);
		}
	}
	for (size_t i = 0; i <= SPLITTING; ++i)
	{
		for (size_t j = 0; j <= SPLITTING; ++j)
		{
			for (size_t p = 0; p < SPLITTING; ++p)
			{
				for (size_t s = 0; s < SPLITTING; ++s)
				{
					//1 треугольник
					overline_u[i][j] -= basicArrays.overline_aa_1[i][j][p][s] * xi[p][s] * u[p][s];
					overline_u[i][j] -= basicArrays.overline_ab_1[i][j][p][s] * xi[p][s] * u[p][s + 1];
					overline_u[i][j] -= basicArrays.overline_ab_1[i][j][p][s] * xi[p][s + 1] * u[p][s];
					overline_u[i][j] -= basicArrays.overline_ac_1[i][j][p][s] * xi[p][s] * u[p + 1][s];
					overline_u[i][j] -= basicArrays.overline_ac_1[i][j][p][s] * xi[p + 1][s] * u[p][s];
					overline_u[i][j] -= basicArrays.overline_bb_1[i][j][p][s] * xi[p][s + 1] * u[p][s + 1];
					overline_u[i][j] -= basicArrays.overline_bc_1[i][j][p][s] * xi[p][s + 1] * u[p + 1][s];
					overline_u[i][j] -= basicArrays.overline_bc_1[i][j][p][s] * xi[p + 1][s] * u[p][s + 1];
					overline_u[i][j] -= basicArrays.overline_cc_1[i][j][p][s] * xi[p + 1][s] * u[p + 1][s];
					//2 треугольник
					overline_u[i][j] -= basicArrays.overline_aa_2[i][j][p][s] * xi[p + 1][s + 1] * u[p + 1][s + 1];
					overline_u[i][j] -= basicArrays.overline_ab_2[i][j][p][s] * xi[p + 1][s + 1] * u[p][s + 1];
					overline_u[i][j] -= basicArrays.overline_ab_2[i][j][p][s] * xi[p][s + 1] * u[p + 1][s + 1];
					overline_u[i][j] -= basicArrays.overline_ac_2[i][j][p][s] * xi[p + 1][s + 1] * u[p + 1][s];
					overline_u[i][j] -= basicArrays.overline_ac_2[i][j][p][s] * xi[p + 1][s] * u[p + 1][s + 1];
					overline_u[i][j] -= basicArrays.overline_bb_2[i][j][p][s] * xi[p][s + 1] * u[p][s + 1];
					overline_u[i][j] -= basicArrays.overline_bc_2[i][j][p][s] * xi[p][s + 1] * u[p + 1][s];
					overline_u[i][j] -= basicArrays.overline_bc_2[i][j][p][s] * xi[p + 1][s] * u[p][s + 1];
					overline_u[i][j] -= basicArrays.overline_cc_2[i][j][p][s] * xi[p + 1][s] * u[p + 1][s];
				}
			}
		}
	}
}
