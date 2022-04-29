/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/

#include "myNM.h"


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

	for (k = 0 ; k < _A.rows; k ++) {
		for (i = k + 1; i < _A.rows; i++) {
			m = _A.at[i][k] / _A.at[k][k];
			
			
			for (j = k ; j < _A.cols; j++) {
				_A.at[i][j] = _A.at[i][j] - ( m * _A.at[k][j] ) ;
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


// Get Solution from LU
void	solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x)
{
	Matrix vecd = createMat(b.rows, b.cols);

	vecd = solveX(P, b);
	
	fwdSub(L, vecd);
	backSub(U, vecd, x);
}


// Multiply Matrix and Vector
Matrix	solveX(Matrix _Ainv, Matrix vecb)
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


// Find Inverse Matrix (Partial Pivoting is applied)
Matrix	inv(Matrix _A)
{
	Matrix Out = createMat(_A.rows, _A.cols);
	Matrix AI = createMat(_A.rows, _A.cols * 2);

	int i, j, k;
	double m;

	Matrix SP = createMat(_A.rows, 1);
	Matrix temp = createMat(1, _A.cols);


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