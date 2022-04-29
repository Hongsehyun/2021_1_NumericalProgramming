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
Matrix	backSub(Matrix _A, Matrix _b)
{
	Matrix Out = createMat(_b.rows, 1);

	// error check: whether _A is a triangular matrix

	//for (int i = _A.rows; i > 0; i--) {
	//	double temp = 0;
	//	for (int j = i + 1; j <= _A.cols; j++)
	//		temp += _A.at[i - 1][j - 1] * Out.at[j - 1][0];
	//	Out.at[i - 1][0] = (_b.at[i - 1][0] - temp) / _A.at[i - 1][i - 1];
	//}

	for (int i = _A.rows-1; i >= 0; i--) {
		double temp = 0;
		for (int j = i + 1; j < _A.cols; j++)
			temp += _A.at[i][j] * Out.at[j][0];
		Out.at[i][0] = (_b.at[i][0] - temp) / _A.at[i][i];
	}

	return Out;
}
