/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Hong Se Hyun
Created          : 03-06-2021
Modified         : 03-06-2021
Language/ver     : C++ in MSVS2019

Description      : [Tutorial]Integration_student.cpp
-------------------------------------------------------------------------------*/

#include "../../include/myNM.h"
//include "myNM.h"
#include "math.h"
#define		pi		3.14159265358979323846264338327950288
#define		Eu		0			// Euler Method
#define		Em		1			// Euler Modified Method

//#include "../../include/myNM.h"


double func(const double t, const double v) {
	double tau = 1;
	double T   = 1 / tau;
	double f   = 10;
	double Vm  = 1;
	double w = 2 * pi * f;

	double OUT = (-v / tau) + ((1 / tau) * Vm * cos(w*t));
	return OUT;
}


int main(int argc, char* argv[])
{
	// PART 1. Euler Method
	printf("\n**************************************************");
	printf("\n                PART 1. Euler Method              ");
	printf("\n**************************************************\n");
	// Initial Value
	double t0		= 0;
	double tf		= 0.1;
	double interval = 0.001;
	double value_EU[101] = { 0,};
	// Euler Method
	printf("Euler Method :: Data of v\n\n");
	odeEU(func, value_EU, t0, tf, interval);
	


	// PART 2. Euler Modified Method
	printf("\n\n\n\n**************************************************");
	printf("\n           PART 2. Euler Modified Method          ");
	printf("\n**************************************************\n");
	// Initial Value
	double value_EM[101] = { 0, };
	// Euler Modified Method
	printf("Euler Modified Method :: Data of v\n\n");
	odeEM(func, value_EM, t0, tf, interval);



	// Extra. You can choose method which you want to run
	printf("\n\n\n\n**************************************************");
	printf("\n                  --  EXTRA --                    \n\n");
	printf("\     You can choose method which you want to run   \n");
	printf("\n**************************************************\n");
	double y[101] = { 0, };
	ode(func, y, t0, tf, interval, 1);			// 0 = Euler Method     1 = Euler Modified Method



	system("pause");
	return 0;
}