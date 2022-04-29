/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Hong Se Hyun
Created          : 03-06-2021
Modified         : 25-06-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H
#define		 PI				3.14159265368979323846264338327950288412

#include "myMatrix.h"

// Extra Assignment
Matrix	planeFit(Matrix _x, Matrix _y, Matrix _z);



extern  void  odePC(double func(const double t, const double v), double y[], double t0, double tf, double h);

// Assignment 9  :: ODE Part 2
// Problem1: Single equation of 1st order ODE
extern	double odeFunc_rc(const double t, const double v);

// Problem2: Single equation of 2nd order ODE
extern	void odeFunc_mck(const double t, const double Y[], double dYdt[]);

// Single Equation : odeEM, odeEU
extern	void odeEU_new(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
extern	void odeEM_new(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
extern	void odeRK2(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
extern	void odeRK4(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);

// 2nd order Equations : sys2RK2, sys2RK4
extern	void sys2RK2(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);
extern	void sys2RK4(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);

extern	void ODEsolver (double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0, int method);



// Assignment 8  :: ODE Part 1
// Euler Method
extern  void  odeEU(double func(const double t, const double v), double y[], double t0, double tf, double h);

// Euler Modified Method
extern  void  odeEM(double func(const double t, const double v), double y[], double t0, double tf, double h);
//extern  void  odeEM(double func(const double t, const double v), double y[], double t0, double tf, double h, double y0);

// Euler Modified Method
extern  void  ode(double func(const double t, const double v), double y[], double t0, double tf, double h, int method);


// Assignment 7  :: Integration
// Rectangular Method for Discrete Data Input
extern	double	IntegrateRect(double _x[], double _y[], int _m);

// Mid Point Method for Discrete Data Input
extern	double	IntegralMid(double _x[], double _y[], int _m);

// Trapezoidal Method for Discrete Data Input
extern	double	trapz(double _x[], double _y[], int _m);

// Simpson 1/3 Method for Function Input
extern	double	integral(double func(const double _x), double _a, double _b, int _n);

// Simpson 3/8 Method for Function Input
extern	double  integral38(double func(const double _x), double _a, double _b, int _n);


// Assignment 6  :: Differentiation
// Return the dy/dx results for the input data. (truncation error: O(h^2))
extern	Matrix	gradient(Matrix _x, Matrix _y);

extern	void	gradient1D(double x[], double y[], double dydx[], int m);

extern	void	printArr(double x[]);

// Return the dy/dx results for the target equation. (truncation error: O(h^2))
extern	Matrix	gradientFunc(double func(const double x), Matrix xin);

// Newton Raphson [ to pass function as input ]
extern	double newtonRaphsonFunc(double func(const double x), double dfunc(const double x), float x0, float tol);


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