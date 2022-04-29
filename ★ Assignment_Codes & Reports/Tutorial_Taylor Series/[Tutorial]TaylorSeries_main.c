/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Jan Park
Created          : 05-03-2021
Modified         : 05-03-2021
Language/ver     : C in MSVS2019

Description      : myNM_main.c
----------------------------------------------------------------*/

#include "myNM.h"

int main(int argc, char* argv[])
{
	double d = 45;
	double x = d * PI / 180;

	double S_N;

	/*===== Select the function to call =====*/
	S_N = sinTaylor(x);
	//S_N = sindTaylor(d);
	//S_N = sinTaylor2(x);
	//S_N = sinTaylor3(x);
	//S_N = sinTaylor4(x);


	printf("\n\n");
	printf("=======================================\n");
	printf("    sin( %f[rad] ) Calculation   \n", x);
	printf("=======================================\n");
	printf("   -  My     result = %3.12f    \n", S_N);
	printf("   -  Math.h result = %3.12f    \n", sin(x));
	printf("   -  absolute err. = %3.12f    \n", S_N - sin(x));
	printf("=======================================\n");

	system("pause");
	return 0;
}