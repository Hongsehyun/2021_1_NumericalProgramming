/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Hong Se Hyun
Created          : 25-06-2021
Modified         : 25-06-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/


#include "myNM.h"
#include <math.h>


Matrix	planeFit(Matrix _x, Matrix _y, Matrix _z) {
	int mx = _x.rows;
	int my = _y.rows;
	int mz = _z.rows;

	Matrix Hong			= createMat(3, 3);
	Matrix Hong_inverse = createMat(3, 3);
	Matrix Sehyun		= createMat(3, 1);
	Matrix matZ			= createMat(3, 1);

	// 1. Error Checking
	if ((mx != my) || (mx != mz) || (my != mz) || mx < 2)
	{
		printf("***********************************************************   \n");
		printf("***********************************************************   \n");
		printf("ERROR: matrix X and Y must be equal and should be more than 1 \n");
		printf("***********************************************************   \n");
		printf("***********************************************************   \n");
	}

	else
	{
		// 2. initialize the values
		double Sxx = 0;
		double Sxy = 0;
		double Sxz = 0;
		double Syy = 0;
		double Syz = 0;
		double Sx = 0;
		double Sy = 0;
		double Sz = 0;

		// 3. Solve for k = 1 to m
		for (int k = 0; k < mx; k++)
		{
			Sxx += _x.at[k][0] * _x.at[k][0];
			Sxy += _x.at[k][0] * _y.at[k][0];
			Sxz += _x.at[k][0] * _z.at[k][0];
			Syy += _y.at[k][0] * _y.at[k][0];
			Syz += _y.at[k][0] * _z.at[k][0];
			Sx  += _x.at[k][0];
			Sy  += _y.at[k][0];
			Sz  += _z.at[k][0];
		}

		// 4. Solve for a2 a1 and a0
		Hong.at[0][0] = Sxx;
		Hong.at[0][1] = Sxy;
		Hong.at[0][2] = Sx;
		Hong.at[1][0] = Sxy;
		Hong.at[1][1] = Syy;
		Hong.at[1][2] = Sy;
		Hong.at[2][0] = Sx;
		Hong.at[2][1] = Sy;
		Hong.at[2][2] = mx;

		Sehyun.at[0][0] = Sxz;
		Sehyun.at[1][0] = Syz;
		Sehyun.at[2][0] = Sz;

		Hong_inverse = inv(Hong);
		matZ = matXvec(Hong_inverse, Sehyun);
	}
	return matZ;
}


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


/*
Matrix	planeFit(Matrix _x, Matrix _y, Matrix _z, Matrix _xq, Matrix _yq, Matrix _zq)
{
	int mx = _x.rows;
	int my = _y.rows;
	int mz = _z.rows;
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
*/


// Interpolation
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




































// Return the dy/dx results for the input data. (truncation error: O(h^2))
Matrix	gradient(Matrix _x, Matrix _y) {

	int m = _x.rows;
	int n = _y.rows;
	Matrix df = createMat(m, 1);

	//Define constant h    ::  h = xi+1 - xi
	double h = _x.at[1][0] - _x.at[0][0];

	//Error Checking
	if ((m != n) || m < 2)
	{
		printf("***********************************************************   \n");
		printf("   ERROR : matrix must be equal and should be more than 1     \n");
		printf("***********************************************************   \n");
	}

	else
	{
		// 1. 3-point Forward difference		// df = [ -3f(x) + 4f(x+1) - f(x+2) ] / 2h
		df.at[0][0] = (-3 * _y.at[0][0] + 4 * _y.at[1][0] - _y.at[2][0]) / (2 * h);

		// 2. 2-point central difference		// df = f(x+1)-f(x-1) / 2h
		for (int i = 1; i < m - 1; i++)
			df.at[i][0] = (_y.at[i + 1][0] - _y.at[i - 1][0]) / (2 * h);

		// 3. 3-point backward difference		// df = [ f(x-2)-4f(x-1)+3f(x) ] / 2h
		df.at[m - 1][0] = (_y.at[m - 3][0] - 4 * _y.at[m - 2][0] + 3 * _y.at[m - 1][0]) / (2 * h);
	}

	return df;
}


// Classical RK4
void sys2RK4(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init) {

	// number of points
	int N = (tf - t0) / h + 1;

	// statr point
	double ti = t0;

	// Initialization
	double K1[2] = { 0 };   // K1 = [K1_y1, K1_y2]    first slope
	double K2[2] = { 0 };   // K2 = [K2_y1, K2_y2]    second slope
	double K3[2] = { 0 };   // K3 = [K3_y1, K3_y2]    third slope
	double K4[2] = { 0 };   // K4 = [K4_y1, K4_y2]    fourth slope

	double Yin[2] = { 0 };
	double K1_y1 = 0;   // y1 :: dydt = z
	double K1_y2 = 0;   // y2 :: dzdt = zdot
	double K2_y1 = 0;
	double K2_y2 = 0;
	double K3_y1 = 0;
	double K3_y2 = 0;
	double K4_y1 = 0;
	double K4_y2 = 0;


	// Initial condition
	y1[0] = y1_init;   //  y(t)
	y2[0] = y2_init;   //  z(t) = dydt(t)

	for (int i = 0; i < N - 1; i++) {

		// Slope 1 [ K1 ]
		Yin[0] = y1[i];		// z
		Yin[1] = y2[i];		// dzdt		

		odeFunc_sys2(ti, Yin, K1);   // Yin = Input     K1 = Output
		K1_y1 = K1[0];
		K1_y2 = K1[1];


		// Slope 2 [ K2 ]
		Yin[0] = y1[i] + ((0.5) * K1_y1 * h);	     	// z              // Euler Method, get Next state(i+1) value 
		Yin[1] = y2[i] + ((0.5) * K1_y2 * h);	    	// dzdt		      // Euler Method, get Next state(i+1) value

		odeFunc_sys2(ti + (0.5 * h), Yin, K2);            // Yin = Input     K2 = Output
		K2_y1 = K2[0];
		K2_y2 = K2[1];


		// Slope 3 [ K3 ]
		Yin[0] = y1[i] + ((0.5) * K2_y1 * h);		   // z              // Euler Method, get Next state(i+1) value 
		Yin[1] = y2[i] + ((0.5) * K2_y2 * h);		   // dzdt		      // Euler Method, get Next state(i+1) value

		odeFunc_sys2(ti + (0.5 * h), Yin, K3);           // Yin = Input     K3 = Output
		K3_y1 = K3[0];
		K3_y2 = K3[1];


		// Slope 4 [ K4 ]
		Yin[0] = y1[i] + (K3_y1 * h);		          // z               // Euler Method, get Next state(i+1) value 
		Yin[1] = y2[i] + (K3_y2 * h);		          // dzdt		      // Euler Method, get Next state(i+1) value

		odeFunc_sys2(ti + h, Yin, K4);                // Yin = Input     K4 = Output
		K4_y1 = K4[0];
		K4_y2 = K4[1];

		// Update
		y1[i + 1] = y1[i] + ((K1_y1 + (2 * K2_y1) + (2 * K3_y1) + K4_y1) * (h / 6));
		y2[i + 1] = y2[i] + ((K1_y2 + (2 * K2_y2) + (2 * K3_y2) + K4_y2) * (h / 6));

		ti += h;
	}
}


// Simpson 1/3 Method for Function Input
double   integral(double func(const double _x), double _a, double _b, int _n) {

	double Out = 0;
	double temp1 = 0;
	double temp2 = 0;

	double h = (_b - _a) / _n;

	for (int i = 1; i < _n; i += 2) {
		double xi = _a + (h * i);			// Interval에 맞춘 x좌표
		temp1 += (4 * func(xi));
	}

	for (int j = 2; j < _n - 1; j += 2) {
		double xi = _a + (h * j);			// Interval에 맞춘 x좌표
		temp2 += (2 * func(xi));
	}

	Out = (h / 3) * (func(_a) + temp1 + temp2 + func(_b));
	return    Out;
}


// Final Exam _ Prob 4.
// Euler Modified Method
void  odeEM(double func(const double t, const double v), double y[], double t0, double tf, double h) {

	// Variable Initializing
	double value_EU = 3;
	double value_EM = 3;
	double slope_1 ;
	double slope_2 ;

	// Print Initial State, Initial Value
	int k = 0;
	y[k] = value_EM;
	//printf("%15.6f\n", y[k]);

	for (double i = 0; i <= 4; i += 1)
	{
		slope_1 = func(i, value_EM);
		value_EU = value_EM + slope_1 * h;

		slope_2 = func(i + h, value_EU);
		value_EM = value_EM + ((slope_1 + slope_2) * (h / 2));

		k++;
		y[k] = value_EM;
		printf("%15.6f\n", y[k]);
	}
}


/*
void  odePC(double func(const double t, const double y), double y[], double t0, double tf, double h) {
	// Variable Initializing
	double value_EU = 3;
	double value_EM = 3;
	double slope_1;
	double slope_2;

	double _tol = 0.00001;
	double ep = 0;

	// Print Initial State, Initial Value
	int k = 0;
	y[k] = value_EM;

	do {

		double ep = y;

		slope_1 = func(k, value_EM);
		value_EU = value_EM + slope_1 * h;

		slope_2 = func(k + h, value_EU);
		value_EM = value_EM + ((slope_1 + slope_2) * (h / 2));

		k++;
		y[k] = value_EM;
		printf("%15.6f\n", y[k]);
	} while (ep > _tol);
}

*/


// Assignment 9  :: ODE Part 2
// Problem1: Single equation of 1st order ODE
// Gradient function for ODE - 1st order single eq.
double odeFunc_rc(const double t, const double v) {
	double tau = 1;
	double T = 1 / tau;
	double f = 10;
	double Vm = 1;
	double w = 2 * PI * f;

	double OUT = (-v / tau) + ((1 / tau) * Vm * cos(w * t));
	return OUT;
}


// Problem2: Single equation of 2nd order ODE
void odeFunc_mck(const double t, const double Y[], double dYdt[])
{
	double m = 1;
	double c = 7;
	double k = 6.9;
	double f = 5;
	double A = 2;

	double Fin = A * cos(2 * PI * f * t);

	// Initial Condition : y(0) = 0.0 m, dy / dt | t = 0 = 0.2 m / s
	// Y[0] == y(t)
	// Y[1] == z(t)
	// dYdt[0] == z(t) == Y[1]
	// dYdt[1] == 2nd order ODE of y :: zdot = (-k*Y - c*Z + Fin)/m;
	dYdt[0] = Y[1];
	dYdt[1] = (-k*Y[0] - c *Y[1] + Fin) / m ;
}


// Single Equation : odeEM, odeEU
void odeEU_new(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0) {

	// Variable Initializing
	double value = 0;
	double slope = 0;

	// Print Initial State, Initial Value
	int k = 0;
	y[k] = value;
	//printf("%15.6f\n", y[k]);

	for (double i = t0; i <= tf; i += h)
	{
		slope = odeFunc(i, value);
		value = value + slope * h;

		k++;
		y[k] = value;
		//printf("%15.6f\n", y[k]);
	}
}


void odeEM_new(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0) {

	// Variable Initializing
	double value_EU = 0;
	double value_EM = 0;
	double slope_1 = 0;
	double slope_2 = 0;

	// Print Initial State, Initial Value
	int k = 0;
	y[k] = value_EM;
	//printf("%15.6f\n", y[k]);

	for (double i = t0; i <= tf; i += h)
	{
		slope_1 = odeFunc(i, value_EM);
		value_EU = value_EM + slope_1 * h;

		slope_2 = odeFunc(i + h, value_EU);
		value_EM = value_EM + ((slope_1 + slope_2) * (h / 2));

		k++;
		y[k] = value_EM;
		//printf("%15.6f\n", y[k]);
	}
}


// General Form of Runge - kutta 2 method uses two slopes to find next value of y.
// General Form of Runge - kutta 2 method set C1 = C2 = 0.5 and alpha = beta = 1. Then It is equal to Euler Modified method.
void odeRK2(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0) {

	double C1 = 0.5;
	double C2 = 0.5;
	double alpha = 1;
	double beta = alpha;  // alpha=beta

	// Numner of points
	int N = ( (tf - t0) / h ) + 1;
	
	// Initialization
	double y_EU = 0;
	double K1 = 0;
	double K2 = 0;
	y[0] = y0;

	for (int i = t0; i < N - 1; i++)
	{
		// First slope
		K1 = odeFunc(t0, y[i]);

		// Second slope
		y_EU = y[i] + ( beta * K1 * h ) ;
		K2 = odeFunc(t0 + alpha * h, y_EU) ;
	
		// Update
		y[i + 1] = y[i] + (C1 * K1 + C2 * K2) * h;
		t0 += h;
	}
}


// General Form of Runge - kutta 4 method uses 4 slopes to find next value of y.
// General Form of Runge - kutta 4 method set C1 = C4 = 1/6 , C2 = C3 = 2 / 6 and alpha = beta = 1/2. 
void odeRK4(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0) {

	double a = 0.5;
	double b = 1;

	// Numner of points
	int N = ((tf - t0) / h) + 1;

	// Initialization
	double y_EU = 0;
	double K1 = 0;
	double K2 = 0;
	double K3 = 0;
	double K4 = 0;
	
	y[0] = y0;

	for (int i = t0; i < N - 1; i++)
	{
		// First slope
		K1 = odeFunc(t0, y[i]);

		// Second slope
		y_EU = y[i] + (a * K1 * h);
		K2 = odeFunc(t0 + (a * h), y_EU);     // In this line, used Euler Mid Method
		
		// Third slope
		y_EU = y[i] + (a * K2 * h);
		K3 = odeFunc(t0 + (a * h), y_EU);     // In this line, used Euler Mid Method

		// Fourth slope
		y_EU = y[i] + K3 * h;
		K4 = odeFunc(t0 + h, y_EU);     // In this line, used Euler Method
		
		// Update
		y[i + 1] = y[i] + (K1 + 2*K2 + 2*K3 + K4) * (h/6) ;
		t0 += h;
	}
}


// ODE RK2:  one of 2nd order ODE <--> two of 1st order ODE
// 2nd order Equations : sys2RK2, sys2RK4
void sys2RK2(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init) {

	// number of points
	int N = (tf - t0) / h + 1;
	
	// statr point
	double ti = t0;

	// Initialization
	double K1[2] = { 0 };   // K1 = [K1_y1, K1_y2]    first slope
	double K2[2] = { 0 };   // K2 = [K2_y1, K2_y2]    second slope
	
	double Yin[2] = { 0 };
	double K1_y1 = 0;   // y1 :: dydt = z
	double K1_y2 = 0;   // y2 :: dzdt = zdot
	double K2_y1 = 0;
	double K2_y2 = 0;

	// Initial condition
	y1[0] = y1_init;   //  y(t)
	y2[0] = y2_init;   //  z(t) = dydt(t)

	for (int i = 0; i < N - 1; i++) {

		// Slope 1 [ K1 ]
		Yin[0] = y1[i];		// z
		Yin[1] = y2[i];		// dzdt		
		
		odeFunc_sys2(ti, Yin, K1);   // Yin = Input     K1 = Output
		K1_y1 = K1[0];
		K1_y2 = K1[1];


		// Slope 2 [ K2 ]
		Yin[0] = y1[i] + K1_y1 * h;		// z              // Euler Method, get Next state(i+1) value 
		Yin[1] = y2[i] + K1_y2 * h;		// dzdt		      // Euler Method, get Next state(i+1) value

		odeFunc_sys2(ti + h, Yin, K2);   // Yin = Input     K2 = Output
		K2_y1 = K2[0];
		K2_y2 = K2[1];


		// Update
		y1[i + 1] = y1[i] + (K1_y1 + K2_y1) * h / 2 ;
		y2[i + 1] = y2[i] + (K1_y2 + K2_y2) * h / 2 ;
		ti += h;
	}
}


void ODEsolver(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0, int method) {

	if (method == 0)
		odeEU_new(odeFunc, y, t0, tf, h, y0);

	else if (method == 1)
		odeEU_new(odeFunc, y, t0, tf, h, y0);

	else if (method == 2)
		odeRK2(odeFunc, y, t0, tf, h, y0);

	else
		odeRK4(odeFunc, y, t0, tf, h, y0);
	
}





/*
// 교수님 코드
void odeEU(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
{
	int N = (tf - t0) / h + 1;
	double ti = t0;
	double slope;

	y[0] = y0;

	for (int i = 0; i < N - 1; i++) {
		slope = odeFunc(ti, y[i]);
		y[i + 1] = y[i] + slope * h;
		ti += h;
	}
}






// Assignment 8  :: ODE Part 1
// Euler Method
void  odeEU(double func(const double t, const double v), double y[], double t0, double tf, double h) {

	// Variable Initializing
	double value = 0;  
	double slope = 0;
	
	// Print Initial State, Initial Value
	int k = 0;
	y[k] = value;
	printf("%15.6f\n", y[k]);
	
	for (double i = t0; i <= tf; i += h)
	{	
		slope = func(i, value);
		value = value + slope * h;

		k++;
		y[k] = value;
		printf("%15.6f\n", y[k]);
	}
}


void  ode(double func(const double t, const double v), double y[], double t0, double tf, double h, int method) {

	if (method == 0)
	{
		odeEU(func, y, t0, tf, h);

	}
	
	if (method == 1)
	{
		odeEM(func, y, t0, tf, h);
	}
}


/*
	for (float i = t0; i <= tf; i += h)
	{

		for (float j = 0; j <= point; j++)
		{
			y[j + 1] = y[j] + (func(i, y[j]) * h);

		}
		//printf("%15.6f\n", y[i]);

	}
*/


// Assignment 7  :: Integration
// Rectangular Method for Discrete Data Input
double IntegrateRect(double _x[], double _y[], int _m) {
	// # of data		: m
	// # of interval	: m-1

	double Out	= 0;
	int		N	= _m - 1;				// Number of Intervals

	for (int i = 0; i < N; i++)
		Out += _y[i] * (_x[i + 1] - _x[i]);

	return Out;
}


// Mid Point Method for Discrete Data Input
double IntegralMid(double _x[], double _y[], int _m) {
	// # of data		: m
	// # of interval	: m-1

	double Out = 0;
	int		N = _m - 1;				// Number of Intervals

	for (int i = 0; i < N; i++)
		Out += ( (_y[i] + _y[i + 1]) / 2 ) * (_x[i + 1] - _x[i]);
	
	return Out;
}


// Trapezoidal Method for Discrete Data Input
double	trapz(double _x[], double _y[], int _m){
	// # of data		: m
	// # of interval	: m-1

	double Out	= 0;
	int		N	= _m - 1;				// Number of Intervals

	for (int i = 0; i < N ; i++)
		Out += ( _y[i] + _y[i + 1] ) * ( _x[i+1] - _x[i] ) ;

	Out = 0.5 * Out;
	return 	Out;
}


// Simpson 3/8 Method for Function Input
double   integral38(double func(const double _x), double _a, double _b, int _n) {

	double Out = 0;
	double temp1 = 0;
	double temp2 = 0;

	double h = (_b - _a) / _n;

	for (int i = 1; i < _n-1; i += 3) {
		double xi = _a + (h * i);			// Interval에 맞춘 x좌표
		temp1 += 3 * ( func(xi) + func(xi+h) );
	}

	for (int j = 3; j < _n - 2; j += 3) {
		double xi = _a + (h * j);			// Interval에 맞춘 x좌표
		temp2 += 2 * func(xi);
	}

	Out = ( (3 * h) / 8 ) * (func(_a) + temp1 + temp2 + func(_b));
	return    Out;
}


/*
// Integration Test Part 2 _ Using Simpson 1/3 Method for Function Input   :: 교수님 풀이
double   integral(double func(const double _x), double _a, double _b, int _n) {
	
	// Out = h/3 * [ y(0) + 4y(i) + 2y(k) + y(n) ]					where i = 1,3,5,7 ... , n-3,n-1		 k = 2,4,6, ... ,n-4,n-2 
	// Out = h/3 * [ y(0) + 4y(i) + 2y(k) + y(n) + 4y(n-1)]			where i = 1,3,5,7 ... , n-3		     k = 2,4,6, ... ,n-4,n-2 
	// Out = h/3 * [ y(0) + 4y(i) + 2y(i+1) + y(n) + 4y(n-1)]		where i = 1,3,5,7 ... , n-3
	
	// i=1 이면 k=2 이고,   i=3 이면 k=4 이므로   i와 k의 관계가 k = i + 1 이라고 할 수 있다.
	// 그리고 i의 마지막 (n-1) 번째 항은 k의 index와 대응되지 않으므로 따로 꺼내줌으로써 i와 k를 하나의 변수에 대해 표현할 수 있도록 해주었다.


	double Out	 = 0;
	double temp  = 0;
	
	double h	 = (_b - _a) / _n;
	
	for (int i = 1 ; i < _n-2 ; i+=2 )
	{
		double xi = _a + ( h * i ) ;			// Interval에 맞춘 x좌표
		temp += ( 4*func(xi) ) + ( 2*func(xi+h) ) ;
	}

	Out = (h / 3) * ( temp + func(_b) + (4*func(_b-h)) );
	return    Out;
}
*/


// Differentiation using Array
void	gradient1D(double x[], double y[], double dydx[], int m) {

	//Define constant h    ::  h = xi+1 - xi
	double h = x[1] - x[0];

	// 1. 3-point Forward difference		// df = [ -3f(x) + 4f(x+1) - f(x+2) ] / 2h
	dydx[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h);

	// 2. 2-point central difference		// df = f(x+1)-f(x-1) / 2h
	for (int i = 1; i < m - 1; i++)
		dydx[i] = (y[i + 1] - y[i - 1]) / (2 * h);

	// 3. 3-point backward difference		// df = [ f(x-2)-4f(x-1)+3f(x) ] / 2h
	dydx[m - 1] = (y[m - 3] - 4 * y[m - 2] + 3 * y[m - 1]) / (2 * h);
}


void	printArr(double x[])
{
	int m = 21;
	for (int i = 0; i < m; i++)
	{
		printf("%15.6f\n", x[i]);
	}
	printf("\n");
}


// Return the dy/dx results for the target equation. (truncation error: O(h^2))
Matrix	gradientFunc(double myFunc(const double x), Matrix xin) {

	Matrix y = createMat(xin.rows , 1);
	Matrix df = createMat(xin.rows , 1);

	// define y[0] to y[n-1]
	for (int i = 0; i < xin.rows; i++)
	{
		y.at[i][0] = myFunc(xin.at[i][0]);
	}

	// define Differentiation of given function
	// Use Gradient() Numerical differentiation
	df = gradient(xin, y);

	return df;
}


// newtonRaphson Code [ to pass function as input ] 
double newtonRaphsonFunc(double func(const double x), double dfunc(const double x), float x0, float tol)
{
	//Assume that f(x) is continuous function.
	//Variable Initializing
	int k = 0;
	int Nmax = 1000;
	double _x = x0;
	double ep = 1000;

	// if f'(_x) is 0, newton Raphson method cannot running.
	// Therefore, we should check if(f'(_x) == 0) 
	if (dfunc(_x) == 0)
	{
		printf("새로운 초기값 x를 입력하세요. \n");
		goto EXIT;     // Immediately close the function newtonRaphson
	}

	//Repeat until(k < Nmax&& ep > _tol)
	do {
		_x = _x - (func(_x) / dfunc(_x));
		ep = fabs(func(_x));

		printf("Iteration:%d \t", k);
		printf("X(n): %f \t", _x);
		printf("Tolerance: %.10f\n", ep);

		k++;

	} while (k < Nmax && ep > tol);

EXIT:
	return _x;
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