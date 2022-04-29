/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H

#include "myMatrix.h"

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
extern	Matrix	solveX(Matrix Ainv, Matrix vecb);

#endif