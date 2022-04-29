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

// Assignment 5 Code
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




// Assignment 2~4
// find F for Jacob
Matrix funcMatrix(double _x, double _y)
{
	Matrix func = createMat(2, 1);

	func.at[0][0] = _y - (pow(2.71828, _x / 2) + pow(2.71828, (-_x) / 2)) / 2;
	func.at[1][0] = 9 * pow(_x, 2) + 25 * pow(_y, 2) - 225;

	//func.at[0][0] = pow(_x, 2) + pow(_y, 2) - 20;
	//func.at[1][0] = _x - _y + 2;

	return func;
}

// find J for Jacob
Matrix jacobMatrix(double _x, double _y)
{
	Matrix dfunc = createMat(2, 2);

	dfunc.at[0][0] = -((pow(2.71828, _x / 2) - pow(2.71828, (-_x) / 2)) / 4);
	dfunc.at[0][1] = 1;
	dfunc.at[1][0] = 18 * _x;
	dfunc.at[1][1] = 50 * _y;

	//dfunc.at[0][0] = _x * 2;
	//dfunc.at[0][1] = _y * 2;
	//dfunc.at[1][0] = 1;
	//dfunc.at[1][1] = -1;

	return dfunc;
}


// newtonRaphson-Jacobian Code
Matrix newtonRaphson_Jacob(Matrix X, Matrix F, Matrix J, double _tol)
{
	Matrix H = createMat(2, 1);
	Matrix L1 = createMat(2, 2);
	Matrix U1 = createMat(2, 2);
	Matrix P1 = createMat(2, 2);

	int k = 1;
	int Nmax = 10;
	double ep = 1000;

	//Repeat until(k < Nmax&& ep > _tol)
	do {
		//ep = 0.1
		F = funcMatrix(X.at[0][0], X.at[1][0]);
		J = jacobMatrix(X.at[0][0], X.at[1][0]);

		for (int i = 0; i < F.rows; i++)
		{
			for (int j = 0; j < F.cols; j++)
			{
				F.at[i][j] = -F.at[i][j];
			}
		}

		//printMat(F, "-F");

		LUdecomp(J, L1, U1, P1);
		solveLU(L1, U1, P1, F, H);

		for (int i = 0; i < X.rows; i++)
		{
			for (int j = 0; j < X.cols; j++)
			{
				X.at[i][j] = X.at[i][j] + H.at[i][j];
			}
		}

		printf("Iteration:%d \n", k);
		printMat(X, "Resulf of Jacob");
		//printf("Tolerance: %.10f\n", ep);
		k++;

	} while (k<Nmax && ep>_tol);

	return X;
}

// Matrix addition
Matrix	addMat(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] + _B.at[i][j];

	return Out;
}


// Forward Substitute
void	fwdSub(Matrix _L, Matrix _d)
{
	int i, j, k;
	double m;

	for (k = 0; k < _L.rows; k++) {
		for (i = k + 1; i < _L.rows; i++) {
			m = _L.at[i][k] / _L.at[k][k];

			for (j = k; j < _L.cols; j++) {
				_L.at[i][j] = _L.at[i][j] - (m * _L.at[k][j]);
			}
			_d.at[i][0] = _d.at[i][0] - m * _d.at[k][0];
		}
	}
}


// Backward Substitute
void	backSub(Matrix _U, Matrix _d, Matrix _x)
{
	for (int i = _U.rows - 1; i >= 0; i--) {
		double temp = 0;
		for (int j = i + 1; j < _U.cols; j++)
			temp += _U.at[i][j] * _x.at[j][0];
		_x.at[i][0] = (_d.at[i][0] - temp) / _U.at[i][i];
	}
}


// Apply Gauss-Elimination without Pivoting
void	gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d)
{
	if (_A.rows != _A.cols) {
		printf("\n*********************************************");
		printf("\n    ! ERROR ! It is not a square Matrix    ");
		printf("\n*********************************************\n");
		printf("\n \n \n \n \n");

		goto EXIT;
	}

	else if (_b.cols != 1) {
		printf("\n*********************************************");
		printf("\n! ERROR ! input Matrix _b is not a column Vector ");
		printf("\n*********************************************\n");
		printf("\n \n \n \n \n");

		goto EXIT;
	}

	else {
		for (int i = 0; i < _A.rows; i++) {
			if (_A.at[i][i] == 0) {
				printf("\n*********************************************");
				printf("\n! ERROR !     Diagonal Components are zero    ");
				printf("\n*********************************************\n");
				printf("\n \n \n \n \n");

				goto EXIT;
			}
		}
	}

	int i, j, k;
	double m;

	for (k = 0; k < _A.rows; k++) {
		for (i = k + 1; i < _A.rows; i++) {
			m = _A.at[i][k] / _A.at[k][k];


			for (j = k; j < _A.cols; j++) {
				_A.at[i][j] = _A.at[i][j] - (m * _A.at[k][j]);
			}
			_b.at[i][0] = _b.at[i][0] - m * _b.at[k][0];
		}
	}

	for (i = 0; i < _A.rows; i++) {
		for (j = 0; j < _A.cols; j++) {
			_U.at[i][j] = _A.at[i][j];
			_d.at[i][0] = _b.at[i][0];
		}
	}

EXIT:
	return;
}


//// Apply LU decomposition without partial Pivoting
//void	LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P)
//{
//
//	// ERROR CHECKING 1
//	if (_A.rows != _A.cols) {
//		printf("\n*********************************************");
//		printf("\n    ! ERROR ! It is not a square Matrix    ");
//		printf("\n*********************************************\n");
//		printf("\n \n \n \n \n");
//	}
//
//	// LU Decomposition
//	// Variable & Matrix Initializing
//	int i, j, k;
//	double m;
//
//	Matrix identity = createMat(_A.rows, _A.cols);
//
//	initMat(_L, 0);
//	initMat(_P, 0);
//
//	identity = eye(identity.rows, identity.cols);
//
//	for (i = 0; i < _A.rows; i++) {
//		for (j = 0; j < _A.cols; j++) {
//			_P.at[i][j] = _P.at[i][j] + identity.at[i][j];
//		}
//	}
//
//	// LU Decomposition CODE without Pivoting
//	for (k = 0; k < _A.rows; k++) {
//		for (k = 0; k < _A.rows; k++) {
//			for (i = k + 1; i < _A.rows; i++) {
//				m = _A.at[i][k] / _A.at[k][k];
//				_L.at[i][k] = m;
//				for (j = k; j < _A.cols; j++) {
//					_A.at[i][j] = _A.at[i][j] - (m * _A.at[k][j]);
//				}
//			}
//		}
//
//		// OUTPUT SETTING
//		for (i = 0; i < _A.rows; i++) {
//			for (j = 0; j < _A.cols; j++) {
//				_U.at[i][j] = _A.at[i][j];
//				_L.at[i][j] = _L.at[i][j] + identity.at[i][j];
//			}
//		}
//	}
//}


// Apply LU decomposition with partial Pivoting
void	LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P)
{
	// LU Decomposition
	// Variable & Matrix Initializing
	int i, j, k;
	double m;

	Matrix identity = createMat(_A.rows, _A.cols);
	Matrix SP = createMat(_A.rows, 1);
	Matrix temp = createMat(1, _A.cols);

	initMat(_L, 0);
	initMat(_P, 0);

	// Making P to identity Matrix
	identity = eye(identity.rows, identity.cols);

	for (i = 0; i < _A.rows; i++) {
		for (j = 0; j < _A.cols; j++) {
			_P.at[i][j] = _P.at[i][j] + identity.at[i][j];
		}
	}

	// LU Decomposition CODE with Pivoting
	for (k = 0; k < _A.rows; k++)
	{
		double max_value = 0;
		double max_row = 0;
		int num = 0;

		// Find Pivot Row (K pivot)
		for (i = k; i < _A.rows; i++) {
			for (j = k; j < _A.cols; j++) {
				if (max_value <= fabs(_A.at[i][j]))
					max_value = fabs(_A.at[i][j]);
			}
			SP.at[i][0] = fabs(_A.at[i][k] / max_value);
			max_value = 0;

			if (max_row < SP.at[i][0])
			{
				max_row = SP.at[i][0];
				num = i;
			}
			initMat(SP, 0);
		}

		// Row Exchange when K and K*pivot is not equal
		if (k != num)
		{
			for (j = 0; j < _A.cols; j++)
			{
				temp.at[0][j] = _P.at[k][j];
				_P.at[k][j] = _P.at[num][j];
				_P.at[num][j] = temp.at[0][j];
			}

			for (j = 0; j < _A.cols; j++)
			{
				temp.at[0][j] = _A.at[k][j];
				_A.at[k][j] = _A.at[num][j];
				_A.at[num][j] = temp.at[0][j];
			}

			for (j = 0; j < _A.cols; j++)
			{
				temp.at[0][j] = _L.at[k][j];
				_L.at[k][j] = _L.at[num][j];
				_L.at[num][j] = temp.at[0][j];
			}
		}

		// Apply Row Reduction
		for (i = k + 1; i < _A.rows; i++) {
			m = _A.at[i][k] / _A.at[k][k];
			_L.at[i][k] = m;
			for (j = k; j < _A.cols; j++) {
				_A.at[i][j] = _A.at[i][j] - (m * _A.at[k][j]);
			}
		}
	}

	// OUTPUT SETTING
	for (i = 0; i < _A.rows; i++) {
		for (j = 0; j < _A.cols; j++) {
			_U.at[i][j] = _A.at[i][j];
			_L.at[i][j] = _L.at[i][j] + identity.at[i][j];
		}
	}
}


// Get Solution from LU
void	solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x)
{
	Matrix vecd = createMat(b.rows, b.cols);

	vecd = matXvec(P, b);

	fwdSub(L, vecd);
	backSub(U, vecd, x);
}


// Find Inverse Matrix (Partial Pivoting is applied)
Matrix	inv(Matrix _A)
{
	Matrix Out = createMat(_A.rows, _A.cols);
	Matrix AI = createMat(_A.rows, _A.cols * 2);
	Matrix SP = createMat(_A.rows, 1);
	Matrix temp = createMat(1, _A.cols);

	int i, j, k;

	// Making [A I] Matrix
	for (i = 0; i < _A.rows; i++)
	{
		for (j = 0; j < _A.cols; j++)
		{
			AI.at[i][j] = _A.at[i][j];

			if (i == j)
			{
				AI.at[i][j + _A.rows] = 1;
			}
			else
			{
				AI.at[i][j + _A.rows] = 0;
			}
		}
	}

	// Gauss Elim with Pivoting
	for (k = 0; k < _A.rows; k++)
	{
		double max_value = 0;
		double max_row = 0;
		int num = 0;

		for (i = k; i < _A.rows; i++) {
			for (j = k; j < _A.cols; j++) {
				if (max_value <= fabs(_A.at[i][j]))
					max_value = fabs(_A.at[i][j]);
			}
			SP.at[i][0] = fabs(_A.at[i][k] / max_value);
			max_value = 0;

			if (max_row < SP.at[i][0])
			{
				max_row = SP.at[i][0];
				num = i;
			}
			initMat(SP, 0);
		}

		if (k != num)
		{
			for (j = 0; j < _A.cols * 2; j++)
			{
				temp.at[0][j] = AI.at[k][j];
				AI.at[k][j] = AI.at[num][j];
				AI.at[num][j] = temp.at[0][j];
			}
		}

		double m1;

		for (i = k + 1; i < _A.rows; i++)
		{
			m1 = AI.at[i][k] / AI.at[k][k];

			for (j = k; j < _A.cols * 2; j++)
			{
				AI.at[i][j] = AI.at[i][j] - (m1 * AI.at[k][j]);
			}
		}
	}

	// Back Substitution
	for (k = _A.rows - 1; k > 0; k--)
	{
		double m2;
		for (i = k - 1; i >= 0; i--)
		{
			m2 = AI.at[i][k] / AI.at[k][k];

			for (j = _A.cols * 2 - 1; j >= k; j--)
			{
				AI.at[i][j] = AI.at[i][j] - (m2 * AI.at[k][j]);
			}
		}
	}

	// Make Diagonal to 1
	for (i = 0; i < _A.rows; i++)
	{
		for (j = _A.cols; j < _A.cols * 2; j++)
		{
			AI.at[i][j] = AI.at[i][j] / AI.at[i][i];
		}
	}

	// Find inverse Matrix of A
	for (i = 0; i < _A.rows; i++) {
		for (j = 0; j < _A.cols; j++) {
			Out.at[i][j] = AI.at[i][j + _A.cols];
		}
	}

	return Out;
}


// Multiply Matrix and Vector
Matrix	matXvec(Matrix _Ainv, Matrix vecb)
{
	Matrix Out = createMat(_Ainv.rows, vecb.cols);

	for (int i = 0; i < _Ainv.rows; i++)
	{
		double temp = 0;
		for (int j = 0; j < _Ainv.cols; j++)
		{
			temp = temp + (_Ainv.at[i][j] * vecb.at[j][0]);
		}
		Out.at[i][0] = temp;
	}

	return Out;
}


// Gauss Jordan을 사용하기 위해서는 Main에 [A b]꼴의 행렬을 받을 수 있는 행렬을 하나 초기화시켜주기!
Matrix	gaussJordan(Matrix _A, Matrix _b)
{
	//Matrix Out = createMat(_A.rows, _A.cols);
	Matrix Ab = createMat(_A.rows, _A.cols + _b.cols);
	Matrix SP = createMat(_A.rows, 1);
	Matrix temp = createMat(1, _A.cols);

	int i, j, k;

	// Making [A b] Matrix
	for (i = 0; i < _A.rows; i++)
	{
		for (j = 0; j < _A.cols; j++)
		{
			Ab.at[i][j] = _A.at[i][j];
			Ab.at[i][_A.cols] = _b.at[i][0];
		}
	}

	// Gauss Elim with Pivoting
	for (k = 0; k < _A.rows; k++)
	{
		double max_value = 0;
		double max_row = 0;
		int num = 0;

		for (i = k; i < _A.rows; i++) {
			for (j = k; j < _A.cols; j++) {
				if (max_value <= fabs(_A.at[i][j]))
					max_value = fabs(_A.at[i][j]);
			}
			SP.at[i][0] = fabs(_A.at[i][k] / max_value);
			max_value = 0;

			if (max_row < SP.at[i][0])
			{
				max_row = SP.at[i][0];
				num = i;
			}
			initMat(SP, 0);
		}

		if (k != num)
		{
			for (j = 0; j < _A.cols + _b.cols; j++)
			{
				temp.at[0][j] = Ab.at[k][j];
				Ab.at[k][j] = Ab.at[num][j];
				Ab.at[num][j] = temp.at[0][j];
			}
		}

		double m1;
		for (i = k + 1; i < _A.rows; i++)
		{
			m1 = Ab.at[i][k] / Ab.at[k][k];

			for (j = k; j <= _A.cols; j++)
			{
				Ab.at[i][j] = Ab.at[i][j] - (m1 * Ab.at[k][j]);
			}
		}
	}

	// Back Substitution
	for (k = _A.rows - 1; k > 0; k--)
	{
		double m2;
		for (i = k - 1; i >= 0; i--)
		{
			m2 = Ab.at[i][k] / Ab.at[k][k];

			for (j = _A.cols; j >= k; j--)
			{
				Ab.at[i][j] = Ab.at[i][j] - (m2 * Ab.at[k][j]);
			}
		}
	}

	// Make Diagonal to 1
	for (i = 0; i < _A.rows; i++)
	{
		Ab.at[i][_A.cols] = Ab.at[i][_A.cols] / Ab.at[i][i];
		for (j = 0; j < _A.cols; j++)
		{
			Ab.at[i][j] = Ab.at[i][j] / Ab.at[i][i];
		}
	}

	return Ab;
}


// Operating QR Decomposition
void QRdecomp(Matrix _A, Matrix _H, Matrix _Q, Matrix _R)
{
	Matrix c = zeros(_A.rows, 1);
	Matrix v = zeros(_A.rows, 1);
	Matrix e = zeros(_A.rows, 1);
	Matrix u = zeros(_A.rows, 1);
	double Norm_of_c = 0;
	double Norm_of_v = 0;
	int i, j, k = 0;

	Matrix memMat_A = zeros(_A.rows, _A.rows);
	Matrix memMat_H = zeros(_A.rows, _A.rows);
	Matrix memMat_Q = zeros(_A.rows, _A.rows);
	Matrix memMat_R = zeros(_A.rows, _A.rows);
	Matrix Transposed_u = zeros(1, _A.rows);

	_R = copyMat(_A);
	_Q = eye(_Q.rows, _Q.cols);

	for (k = 0; k < _R.rows - 1; k++)
	{
		//initialize values for each K
		for (i = 0; i < _H.rows; i++)
		{
			for (j = 0; j < _H.cols; j++)
			{
				if (i == j)
					_H.at[i][j] = 1;
				else
					_H.at[i][j] = 0;
			}

			c.at[i][0] = 0;
			v.at[i][0] = 0;
			e.at[i][0] = 0;
		}


		for (i = k; i < _R.rows; i++)
		{
			c.at[i][0] = _R.at[i][k];
		}

		if (c.at[k][0] >= 0)
		{
			e.at[k][0] = 1;
		}

		else
		{
			e.at[k][0] = -1;
		}


		// make vector v using C and norm C
		// make vector u using v and norm V
		Norm_of_c = Norm(c);
		for (i = 0; i < _R.rows; i++)
		{
			v.at[i][0] = c.at[i][0] + (Norm_of_c * e.at[i][0]);
		}

		Norm_of_v = Norm(v);
		for (i = 0; i < _R.rows; i++)
		{
			u.at[i][0] = v.at[i][0] / Norm_of_v;
		}

		//Transpose u
		//Transposed_u = transpose(u);
		for (i = 0; i < Transposed_u.rows; i++)
		{
			for (j = 0; j < Transposed_u.cols; j++)
			{
				Transposed_u.at[i][j] = u.at[j][i];
			}
		}

		//Find Household and Reflection Matrix H
		multiply_Mat(u, Transposed_u, memMat_H);

		for (i = 0; i < _A.rows; i++)
		{
			for (j = 0; j < _A.rows; j++)
			{
				memMat_H.at[i][j] = 2 * (memMat_H.at[i][j]);
			}
		}

		//we already initailized H as I, just subtract H and memMat !
		for (i = 0; i < _H.rows; i++)
		{
			for (j = 0; j < _H.cols; j++)
			{
				_H.at[i][j] -= memMat_H.at[i][j];
			}
		}

		// Using H, we can find Q, R
		multiply_Mat(_Q, _H, memMat_Q);
		for (i = 0; i < memMat_Q.rows; i++)
		{
			for (j = 0; j < memMat_Q.cols; j++)
			{
				_Q.at[i][j] = memMat_Q.at[i][j];
			}
		}

		multiply_Mat(_H, _R, memMat_R);
		for (i = 0; i < memMat_R.rows; i++)
		{
			for (j = 0; j < memMat_R.cols; j++)
			{
				_R.at[i][j] = memMat_R.at[i][j];
			}
		}

		multiply_Mat(_R, _Q, memMat_A);
		for (i = 0; i < memMat_A.rows; i++)
		{
			for (j = 0; j < memMat_A.cols; j++)
			{
				_A.at[i][j] = memMat_A.at[i][j];
			}
		}
	}
	printMat(_H, "H");
	printMat(_Q, "Q");
	printMat(_R, "R");
	printMat(_A, "A");
}


Matrix	eig(Matrix _A, Matrix _H, Matrix _Q, Matrix _R)
{

	// ERROR CHECKING 1
	if (_A.rows != _A.cols) {
		printf("\n*********************************************");
		printf("\n    ! ERROR ! It is not a square Matrix    ");
		printf("\n*********************************************\n");
		printf("\n \n \n \n \n");
	}

	Matrix eigenValue = zeros(_A.rows, 1);

	for (int z = 0; z <= 100; z++)
	{

		Matrix c = zeros(_A.rows, 1);
		Matrix v = zeros(_A.rows, 1);
		Matrix e = zeros(_A.rows, 1);
		Matrix u = zeros(_A.rows, 1);
		double Norm_of_c = 0;
		double Norm_of_v = 0;
		int i, j, k = 0;

		Matrix memMat_A = zeros(_A.rows, _A.rows);
		Matrix memMat_H = zeros(_A.rows, _A.rows);
		Matrix memMat_Q = zeros(_A.rows, _A.rows);
		Matrix memMat_R = zeros(_A.rows, _A.rows);
		Matrix Transposed_u = zeros(1, _A.rows);

		_R = copyMat(_A);
		_Q = eye(_Q.rows, _Q.cols);

		for (k = 0; k < _R.rows - 1; k++)
		{
			//initialize values for each K
			for (i = 0; i < _H.rows; i++)
			{
				for (j = 0; j < _H.cols; j++)
				{
					if (i == j)
						_H.at[i][j] = 1;
					else
						_H.at[i][j] = 0;
				}

				c.at[i][0] = 0;
				v.at[i][0] = 0;
				e.at[i][0] = 0;
			}


			for (i = k; i < _R.rows; i++) {
				c.at[i][0] = _R.at[i][k];
			}

			if (c.at[k][0] >= 0) {
				e.at[k][0] = 1;
			}

			else {
				e.at[k][0] = -1;
			}

			// make vector v using C and norm C
			// make vector u using v and norm V
			Norm_of_c = Norm(c);
			for (i = 0; i < _R.rows; i++) {
				v.at[i][0] = c.at[i][0] + (Norm_of_c * e.at[i][0]);
			}

			Norm_of_v = Norm(v);
			for (i = 0; i < _R.rows; i++) {
				u.at[i][0] = v.at[i][0] / Norm_of_v;
			}

			//Transpose u
			for (i = 0; i < Transposed_u.rows; i++) {
				for (j = 0; j < Transposed_u.cols; j++) {
					Transposed_u.at[i][j] = u.at[j][i];
				}
			}

			//Find Household and Reflection Matrix H
			multiply_Mat(u, Transposed_u, memMat_H);

			for (i = 0; i < _A.rows; i++) {
				for (j = 0; j < _A.rows; j++) {
					memMat_H.at[i][j] = 2 * (memMat_H.at[i][j]);
				}
			}

			//we already initailized H as I, just subtract H and memMat !
			for (i = 0; i < _H.rows; i++) {
				for (j = 0; j < _H.cols; j++) {
					_H.at[i][j] -= memMat_H.at[i][j];
				}
			}

			// Using H, we can find Q, R
			multiply_Mat(_Q, _H, memMat_Q);
			for (i = 0; i < memMat_Q.rows; i++) {
				for (j = 0; j < memMat_Q.cols; j++) {
					_Q.at[i][j] = memMat_Q.at[i][j];
				}
			}

			multiply_Mat(_H, _R, memMat_R);
			for (i = 0; i < memMat_R.rows; i++) {
				for (j = 0; j < memMat_R.cols; j++) {
					_R.at[i][j] = memMat_R.at[i][j];
				}
			}

			multiply_Mat(_R, _Q, memMat_A);
			for (i = 0; i < memMat_A.rows; i++) {
				for (j = 0; j < memMat_A.cols; j++) {
					_A.at[i][j] = memMat_A.at[i][j];
				}
			}
		}
	}

	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++) {
			if (i == j)
				eigenValue.at[i][0] = _A.at[i][j];
		}
	}

	return eigenValue;
}


double	condNumber(Matrix _A)
{
	Matrix transposeA = createMat(_A.cols, _A.rows);
	Matrix ATA = createMat(_A.cols, _A.cols);
	Matrix tempA = createMat(_A.cols, _A.cols);
	Matrix eigenATA = createMat(ATA.rows, ATA.cols);

	Matrix ATAmatH = createMat(ATA.rows, ATA.cols);
	Matrix ATAmatQ = createMat(ATA.rows, ATA.cols);
	Matrix ATAmatR = createMat(ATA.rows, ATA.cols);

	double maxeig_ATA = 0;
	double mineig_ATA = 0;
	double cond = 0;


	//Find ATA
	transposeA = transpose(_A);
	multiply_Mat(transposeA, _A, tempA);
	ATA = copyMat(tempA);

	//Find EigenValue of A
	eigenATA = eig(ATA, ATAmatH, ATAmatQ, ATAmatR);

	//Find max EigenValue of A
	maxeig_ATA = max(eigenATA);
	mineig_ATA = min(eigenATA);

	// ERROR CHECKING 2
	if (mineig_ATA == 0) {
		printf("\n*********************************************");
		printf("\n     ! ERROR ! caused division by ZERO     ");
		printf("\n*********************************************\n");
		printf("\n \n \n \n \n");
	}

	cond = sqrt(maxeig_ATA / mineig_ATA);

	return cond;
}


/*
Matrix	eig(Matrix _A, Matrix _H, Matrix _Q, Matrix _R)
{
	Matrix UpdateA = zeros(_A.rows, _A.cols);
	Matrix eigenValue = zeros(_A.rows, 1);
	Matrix temp = zeros(_A.rows, _A.cols);

	Matrix memMat_H = zeros(_A.rows, _A.rows);
	Matrix memMat_Q = zeros(_A.rows, _A.rows);
	Matrix memMat_R = zeros(_A.rows, _A.rows);

	for (int i = 0; i < UpdateA.rows; i++) {
		for (int j = 0; j < UpdateA.cols; j++) {
			UpdateA.at[i][j] = _A.at[i][j];
		}
	}

	// error checking 1
	if (_A.rows != _A.cols) {
		printf("\n*********************************************");
		printf("\n    ! error ! it is not a square matrix    ");
		printf("\n*********************************************\n");
		printf("\n \n \n \n \n");
	}


	for (int k = 0; k < 100 ; k++)
	{
		QRdecomp(UpdateA, _H, _Q, _R);

		multiply_Mat(_R, _Q, temp);

		for (int i = 0; i < UpdateA.rows; i++) {
			for (int j = 0; j < UpdateA.cols; j++) {
				UpdateA.at[i][j] = temp.at[i][j];
			}
		}

		initMat(_R, 0);
		initMat(_Q, 0);
		initMat(_H, 0);
	}


	for (int i = 0; i < _A.rows; i++)
	{
		for (int j = 0; j < _A.cols; j++)
		{
			if (i == j)
				eigenValue.at[i][0] = UpdateA.at[i][j];
		}
	}

	return eigenValue;
}
*/