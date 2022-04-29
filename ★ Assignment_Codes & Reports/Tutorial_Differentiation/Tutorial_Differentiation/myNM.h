/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Hong Se Hyun
Created          : 09-05-2021
Modified         : 09-05-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H

#include "myMatrix.h"

// Assignment 5
extern Matrix	linearFit(Matrix _x, Matrix _y);

extern Matrix	linearInterp(Matrix _x, Matrix _y, Matrix _xq);

extern Matrix	arr2Mat(double* _1Darray, int _rows, int _cols);

// Matrix addition
extern	Matrix	addMat(Matrix _A, Matrix _B);

// Forward Substitute
extern	void	fwdSub(Matrix _L, Matrix _d);

// Backward Substitute
extern	void	backSub(Matrix _U, Matrix _d, Matrix _x);

// GaussElim without Pivoting
extern	void	gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d);

// LU decomposition without Pivoting
extern	void	LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P);

// Solve LU decomposition
extern	void	solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x);

// Inverse Matrix
extern	Matrix	inv(Matrix _A);

// Multiply Matrix
extern	Matrix	matXvec(Matrix Ainv, Matrix vecb);

// Gauss-Jordan Elimination
extern	Matrix	gaussJordan(Matrix _A, Matrix _b);

// QR Decomposition
extern	void	QRdecomp(Matrix _A, Matrix _H, Matrix _Q, Matrix _R);

// Find EigenValues from QR Decomposition
// extern	Matrix	eigenvalue(Matrix _A);
extern	Matrix	eig(Matrix _A, Matrix _H, Matrix _Q, Matrix _R);

// Condition Number
extern	double	condNumber(Matrix _A);

//Given Function
extern Matrix funcMatrix(double _x, double _y);
extern Matrix jacobMatrix(double _x, double _y);

//NewtonRhapson function
extern Matrix	newtonRaphson_Jacob(Matrix X, Matrix F, Matrix J, double _tol);

#endif