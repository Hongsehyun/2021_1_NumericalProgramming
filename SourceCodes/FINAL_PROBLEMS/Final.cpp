/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Jan Park, YKKIM
Created          : 2021-06-03
Modified         : 2021-06-03  by YKKIM
Language/ver     : C++ in MSVS2017

Description      : [Tutorial]ODE_IVP_student.c
-------------------------------------------------------------------------------*/

#include "../../include/myNM.h"
//#include "myNM.h"

#include <stdio.h>
#include <math.h>

#define	PI	   3.14159265368979323846264338327950288412
#define ODE_EU 0
#define ODE_EM 1
#define ODE_RK2 2
#define ODE_RK4 3


double cantilever(const double x) {

	double w = 10000;		//[N / m]
	double L = 1;				//[m]
	double E = 200000000000;	//[Pa]
	double I = 0.00021;		// [m ^ 4]
	
	double U = pow((-w * (pow(x, 2) / 2)), 2) / (2 * E * I);
	return U;
}


void ode_pendulum(const double t, const double Y[], double dYdt[]) {

	double c = 0.16;	     // [N * s / m
	double m = 0.5;		 	 // [kg]
	double L = 1.2;			 // [m]
	double g = 9.8;			 // [m / s ^ 2]

	// Y[0] == y(t)
	// Y[1] == z(t)
	// dYdt[0] == z(t) == Y[1]	
	// dYdt[1] == 2nd order ODE of y :: zdot = (-k*Y - c*Z + Fin)/m;
	
	double fin = ((-c / m)*Y[1]) - ((g / L) * sin(Y[0]));

	dYdt[0] = Y[1];
	dYdt[1] = fin;
}



int main(int argc, char* argv[])
{
	// Curve Fitting . Prob 1.
	/*
	double V_array[] = { 0, 2.5, 6.2, 8.2, 9.6, 12.5, 13.1 , 14.1 , 16.6, 18.3 };
	double I_array[] = { 0,2,5,7,8,10,11,12,14,15 };

	Matrix V = createMat(10, 1);
	Matrix I = createMat(10, 1);
	V = arr2Mat(V_array, 10, 1);
	I = arr2Mat(I_array, 10, 1);

	Matrix coefficientZ = createMat(2, 1);
	coefficientZ = linearFit(I, V);

	double a0 = coefficientZ.at[0][0];
	double a1 = coefficientZ.at[1][0];
	//double inputTemp1 = 15;
	//double predict_R = (a0 * inputTemp1) + a1;

	printf("\n-------------------------------------------------------------");
	printf("\n         Coefficient              ");
	printf("\n-------------------------------------------------------------\n");
	printMat(coefficientZ, "Coefficient [a0] and a[1]");

	printf("Equation Form : f(x) = %fx + %f \n\n", a0, a1);
	//printf("At voltate %f, the resistor value is  %f \n\n\n", inputTemp1, predict_R);
	*/


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//


	// Curve Fitting . Prob 2.
	/*
	double U_array[] = { 4.72, 12.49, 20.03, 28.33, 37.47, 41.43, 48.38, 55.06, 66.77, 59.16, 54.45, 47.21, 42.75, 32.71, 25.43, 8.18 };
	double V_array[] = { 7.18, 7.3, 7.37, 7.42, 7.47, 7.5, 7.53, 7.55, 7.58, 7.56 ,7.55, 7.53, 7.51, 7.47, 7.44, 7.28 };
	Matrix U = createMat(16, 1);
	Matrix V = createMat(16, 1);
	Matrix lnU = createMat(16, 1);

	U = arr2Mat(U_array, 16, 1);
	V = arr2Mat(V_array, 16, 1);

	//double mod_y_array[] = { 0, };
	//double lny = 0;

	for (int i = 0; i <= 15; i++) {
		lnU.at[i][0] = log(U.at[i][0]);
	}

	Matrix coefficientZ = createMat(2, 1);
	coefficientZ = linearFit(V, lnU);

	double a0 = coefficientZ.at[0][0];
	double a1 = exp(coefficientZ.at[1][0]);

	printMat(coefficientZ, "Coefficient [a0] and a[1]");
	*/


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//


	// differentiation . Prob 3.
	/*
	double k = 240;
	double X_Array[] = { 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.01 };             // 단위 : m
	double T_Array[] = { 473, 446.3, 422.6, 401.2, 382, 364.3, 348.0, 332.7, 318.1, 304.0, 290.1 };   // 단위 : K

	Matrix flux = createMat(2, 1);
	Matrix dataX = createMat(11, 1);
	Matrix dataT = createMat(11, 1);
	Matrix Result = createMat(11, 1);

	dataX = arr2Mat(X_Array, 11, 1);
	dataT = arr2Mat(T_Array, 11, 1);

	Result = gradient(dataX, dataT);

	flux.at[0][0] = (-k*Result.at[0][0])/1000;           // 단위 : kW/m^2
	flux.at[1][0] = (-k*Result.at[10][0])/1000;

	printMat(Result, "Result");
	printMat(flux, "flux");
	*/


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//


	// Integration . Prob 4
	/*
	double start     = 0;	// 적분구간 a
	double finish	 = 1;	// 적분구간 b
	double n		 = 30;	// Number of Intervals = 30   // if using 38 simpson method, intervals should be multiplication form of 3

	double U = 0;

	U = integral38(cantilever, start, finish, n);
	printf("Resulf of Inteegral is :: %f\n", U);
	*/


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//


	// ODE . Prob 5
	// Parameter Definitions
	double t0 = 0;
	double tf = 9;
	double h = 0.1;
	double N = (tf - t0) / h + 1;

	double d1_sysRK4[200] = { 0 };
	double d2_sysRK4[200] = { 0 };


	// Initial values
	double theta0 = PI / 2;
	double w0 = 0.2;
	
	sys2RK4(ode_pendulum, d1_sysRK4, d2_sysRK4, t0, tf, h, theta0, w0);

	// Print outputs
	printf("/*---------------------------*/\n");
	printf("/ 2nd Order ODE : MCK example /\n");
	printf("/*---------------------------*/\n");
	printf(" - Total number of data N=%d \n", N);
	for (int i = 0; i < N; i++) {
		printf("RK4: t= %f\tRKy= %f\tRKv= %f\n ", t0 + i * h, d1_sysRK4[i], d2_sysRK4[i]);
	}
	



	system("pause");
	return 0;
}