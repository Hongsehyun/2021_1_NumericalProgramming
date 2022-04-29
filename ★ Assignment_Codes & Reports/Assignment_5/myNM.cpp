/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Hong Se Hyun
Created          : 09-05-2021
Modified         : 09-05-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/

#include "myNM.h"
#include <math.h>

// Returns the parameters of the linear least square function.
Matrix	linearFit(Matrix _x, Matrix _y) {

	int mx = _x.rows;
	int my = _y.rows;
	double a1, a0;

	Matrix matZ = createMat(2, 1);

	// 1. Error Checking
	if ((mx != my) || mx < 2)
	{
		printf("***********************************************************   \n");
		printf("ERROR: matrix X and Y must be equal and should be more than 1 \n");
		printf("***********************************************************   \n");
		a1 = 0;
		a0 = 0;
		// goto then exit
	}

	else
	{
		// 2. initialize the values(Sx, Sxx, Sy, Sxy)
		double Sx = 0;
		double Sxx = 0;
		double Sxy = 0;
		double Sy = 0;

		// 3. Solve Sx, Sxx, Sy, Sxy, for k = 1 to m
		for (int k = 0; k < mx; k++)
		{
			Sx += _x.at[k][0];
			Sxx += _x.at[k][0] * _x.at[k][0];
			Sxy += _x.at[k][0] * _y.at[k][0];
			Sy += _y.at[k][0];
		}

		// 4. Solve for a0 and a1
		double DEN = (mx * Sxx - Sx * Sx);
		a1 = (mx * Sxy - Sx * Sy) / DEN;
		a0 = (Sxx * Sy - Sxy * Sx) / DEN;
	}
	// 5. Return z = [a1, a0]
	matZ.at[0][0] = a1;
	matZ.at[1][0] = a0;

	return matZ;
}


Matrix	linearInterp(Matrix _x, Matrix _y, Matrix _xq)
{
	int mx = _x.rows;
	int my = _y.rows;
	double a1 = 0;
	double a0 = 0;

	Matrix mem_coeffi = createMat(2, _x.rows - 1);
	Matrix matZ = createMat(_xq.rows, 1);

	// 1. Error Checking
	if ((mx != my) || mx < 2)
	{
		printf("***********************************************************   \n");
		printf("ERROR: matrix X and Y must be equal and should be more than 1 \n");
		printf("***********************************************************   \n");
		goto EXIT;
	}

	else
	{	
		for (int i = 0; i < _xq.rows; i++)
		{
			for (int j = 0; j < _x.rows - 1; j++)
			{
				// 2. Define Coefficient [a1] and [a0] for each function by using Lagrange Form
				double DEN = _x.at[j][0] - _x.at[j + 1][0];

				a1 = (_y.at[j][0] - _y.at[j + 1][0]) / DEN;
				a0 = (_y.at[j + 1][0] * _x.at[j][0] - _y.at[j][0] * _x.at[j + 1][0]) / DEN;

				mem_coeffi.at[0][j] = a1;
				mem_coeffi.at[1][j] = a0;

				// Condition :: _x.at[j][0] <= _xq.at[i][0] <=_x.at[j+1][0] 
				if (_x.at[j][0] <= _xq.at[i][0] && _xq.at[i][0] <= _x.at[j + 1][0])
				{
					// 3. make each functions :: 'f(x) = a1x + a0' form
					matZ.at[i][0] = mem_coeffi.at[0][j] * _xq.at[i][0] + mem_coeffi.at[1][j];
				}
			}
		}
	}

EXIT:
	return matZ;
}


// Create a matrix from 1D-array
Matrix	arr2Mat(double* _1Darray, int _rows, int _cols)
{
	Matrix Output = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			Output.at[i][j] = _1Darray[i * _cols + j];

	return Output;
}


/*
Matrix	linearInterp(Matrix _x, Matrix _y, Matrix _xq)
{
	int mx = _x.rows;
	int my = _y.rows;
	double a1 = 0;
	double a0 = 0;

	Matrix matZ = createMat(_xq.rows, 1);

	if ((mx != my) || mx < 2)
	{
		printf("***********************************************************   \n");
		printf("ERROR: matrix X and Y must be equal and should be more than 1 \n");
		printf("***********************************************************   \n");
		a1 = 0;
		a0 = 0;
	}

	else
	{
		for (int i = 0; i < mx - 1; i++)
		{
			double DEN = _x.at[i][0] - _x.at[i + 1][0];

			a1 = (_y.at[i][0] - _y.at[i + 1][0]) / DEN;
			a0 = (_y.at[i + 1][0] * _x.at[i][0] - _y.at[i][0] * _x.at[i + 1][0]) / DEN;

			matZ.at[2 * i][0] = a1 * _xq.at[2 * i][0] + a0;
			matZ.at[2 * i + 1][0] = a1 * _xq.at[2 * i + 1][0] + a0;
		}

		matZ.at[(2 * mx) - 2][0] = a1 * _xq.at[(2 * mx) - 2][0] + a0;

		return matZ;
	}
}
*/



