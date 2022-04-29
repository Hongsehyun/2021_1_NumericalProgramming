/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Hong Se Hyun
Created          : 18-05-2021
Modified         : 18-05-2021
Language/ver     : C++ in MSVS2019

Description      : [Tutorial]Integration_student.cpp
-------------------------------------------------------------------------------*/

#include "myNM.h"
#include "math.h"

double myFunc(const double x) {
	double y = sqrt(1 - pow(x, 2));
	return y;
}



int main(int argc, char* argv[])
{
	// PART 1. Integration from Datasets
	printf("\n**************************************************");
	printf("\n        PART 1. Integration from Datasets         ");
	printf("\n**************************************************\n");

	double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
	double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
	int M = sizeof(x) / sizeof(x[0]);

	// Rectangular Method for Discrete Data Input
	double I_rect = IntegrateRect(x, y, M);
	printf("I_rect  = %f\n\n", I_rect);

	// Mid Point Method for Discrete Data Input
	double I_midpoint = IntegralMid(x, y, M);
	printf("I_midpoint = %f\n\n", I_midpoint);

	// Trapezoidal Method for Discrete Data Input
	double I_trapz = trapz(x, y, M);
	printf("I_trapz = %f\n\n", I_trapz);



	// PART 2. Integration from a Function
	printf("\n**************************************************");
	printf("\n        PART 2. Integration from a Function       ");
	printf("\n**************************************************\n");

	double a = -1;		// 利盒备埃 a
	double b = 1;		// 利盒备埃 b
	double n = 12;      // Number of Intervals = 12

	// Simpson 1/3 Method for Function Input
	double I_simpson13 = integral(myFunc, a, b, n);
	printf("I_simpson1/3  = %f\n\n", I_simpson13);

	// Simpson 3/8 Method for Function Input
	double I_simpson38 = integral38(myFunc, a, b, n);
	printf("I_simpson3/8  = %f\n\n", I_simpson38);



	system("pause");
	return 0;
}