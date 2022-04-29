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


// Apply back-substitution
void	backSub(Matrix _U, Matrix _d, Matrix _x)
{
	for (int i = _U.rows - 1; i >= 0; i--) {
		double temp = 0;
		for (int j = i + 1; j < _U.cols; j++)
			temp += _U.at[i][j] * _x.at[j][0];
		_x.at[i][0] = (_d.at[i][0] - temp) / _U.at[i][i];
	}
}


//Apply Gauss-Elimination without Pivoting
void	gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d)
{

	// ERROR CHECKING 1
	if (_A.rows != _A.cols) {
		printf("\n*********************************************");
		printf("\n    ! ERROR ! It is not a square Matrix    ");
		printf("\n*********************************************\n");
		printf("\n \n \n \n \n");

		goto EXIT;		// immediately close the function 'gaussElim'
	}

	// ERROR CHECKING 2
	else if (_b.cols != 1) {
		printf("\n*********************************************");
		printf("\n! ERROR ! input Matrix _b is not a column Vector ");
		printf("\n*********************************************\n");
		printf("\n \n \n \n \n");

		goto EXIT;		// immediately close the function 'gaussElim'
	}

	// ERROR CHECKING 3
	else {
		for (int i = 0; i < _A.rows; i++) {
			if (_A.at[i][i] == 0) {
				printf("\n*********************************************");
				printf("\n! ERROR !     Diagonal Components are zero    ");
				printf("\n*********************************************\n");
				printf("\n \n \n \n \n");

				goto EXIT;		// immediately close the function 'gaussElim'
			}
		}
	}



	// Gauss Elimination
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