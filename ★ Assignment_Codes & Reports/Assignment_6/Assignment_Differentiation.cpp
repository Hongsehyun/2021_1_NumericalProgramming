/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Hong se hyun
Created          : 12-05-2021
Modified         : 12-05-2021
Language/ver     : C++ in MSVS2019

Description      : [Tutorial]Differentiation_student.cpp
-------------------------------------------------------------------------------*/

#define Assignment   6		   // enter your assignment number
#define eval		 0         // set

#include "myNM.h"

// Define a function that defines the target equation.
double myFunc(const double x) {
	double y = x * x * x;
	return y;
}

double mydFunc(const double x) {
	double y = 3 * x * x;
	return y;
}



int main(int argc, char* argv[])
{
	/*    [※ DO NOT EDIT IT !!!]   Resources file path setting for evaluation   */
	std::string path = "C:/NM_data_2021/Assignment" + std::to_string(Assignment) + "/";

#if eval
	path += "eval/";
#endif

	/*==========================================================================*/
	//Matrix t = txt2Mat(path, "Q1_vect");
	//Matrix pos = txt2Mat(path, "Q1_vecx");
	//Matrix xin = txt2Mat(path, "Q2_vecxin");
	/*==========================================================================*/


	/*==========================================================================*/
	// txt2Mat을 사용하지 않고, Array 와 Arr2Mat를 이용하여 Data Setting
	int m = 21;

	double t_arr[21] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6 ,3.8 ,4.0 };
	double pos_arr[21] = { -5.87, -4.23, -2.55, -0.89, 0.67, 2.09, 3.31, 4.31, 5.06, 5.55, 5.78, 5.77, 5.52, 5.08, 4.46, 3.72, 2.88, 2.0, 1.1, 0.23, -0.59 };
	double vel_arr[21] = { 0.0 };
	double acc_arr[21] = { 0.0 };

	Matrix t	= arr2Mat(t_arr, 21, 1);
	Matrix xin	= arr2Mat(t_arr, 21, 1);
	Matrix pos	= arr2Mat(pos_arr, 21, 1);
	/*==========================================================================*/


	// PART 1-1 :: gradient with Matrix
	printf("\n**************************************************");
	printf("\n|                    PART 1.1                    |");
	printf("\n|                                                |");
	printf("\n|            - Gradient with Matrix -            |");
	printf("\n**************************************************\n");

	Matrix vel = gradient(t, pos);
	Matrix acc = gradient(t, vel);

	printMat(t, "t");
	printMat(pos, "pos");
	printMat(vel, "vel");
	printMat(acc, "acc");


	// PART 1-2 :: gradient with Array
	printf("\n**************************************************");
	printf("\n|                     PART 1.2                   |");
	printf("\n|                                                |");
	printf("\n|             - Gradient with Array -            |");
	printf("\n**************************************************\n");
		
	gradient1D(t_arr, pos_arr, vel_arr, m);
	gradient1D(t_arr, vel_arr, acc_arr, m);

	printf("vel_arr = \n");
	printArr(vel_arr);
	printf("acc_arr = \n");
	printArr(acc_arr);

	
	// PART 2-1 :: gradient [ From a Function ]
	printf("\n**************************************************");
	printf("\n|                     PART 2.1                   |");
	printf("\n|                                                |");
	printf("\n|         - Gradient from the Function -         |");
	printf("\n**************************************************\n");

	Matrix dydx = gradientFunc(myFunc, xin);

	printMat(xin, "xin");
	printMat(dydx, "dydx");

	
	// PART 2-2 :: Application to Newton-Raphson
	printf("\n**************************************************");
	printf("\n|                    PART 2.2                   |");
	printf("\n|                                               |");
	printf("\n|          - Newton-Raphson Update Ver. -       |");
	printf("\n**************************************************\n");

	double _tol = 0.0000000000000000001;
	double _x = 2.0;						 // try to choose 'x' close true value
	double NRM_result;

	printf("Newton-Raphson Method:\n");
	NRM_result = newtonRaphsonFunc(myFunc, mydFunc, _x, _tol);

	printf("Final Solution: %f \t", NRM_result);
	printf("\n");	



	system("pause");
	return 0;
}