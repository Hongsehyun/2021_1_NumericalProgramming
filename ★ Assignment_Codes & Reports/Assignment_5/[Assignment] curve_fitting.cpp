/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Hong Se Hyun
Created          : 09-05-2021
Modified         : 09-05-2021
Language/ver     : C++ in MSVS2019

Description      : [Tutorial]curve_fitting.cpp
-------------------------------------------------------------------------------*/

#define Assignment   5		   // enter your assignment number
#define eval		 0         // set

#include "myNM.h"

int main(int argc, char* argv[])
{
	/*    [¡Ø DO NOT EDIT IT !!!]   Resources file path setting for evaluation   */
	std::string path = "C:/NM_data_2021/Assignment" + std::to_string(Assignment) + "/";

#if eval
	path += "eval/";
#endif

	/*==========================================================================*/
	/*               Variables declaration & initialization					    */
	/*--------------------------------------------------------------------------*/
	/*   - You can change the variable names							        */
	/*   - However, you must use the specified file name					    */
	/*      : For each assignment, the file name will be notified on HISNET     */
	/*==========================================================================*/
	Matrix vecT = txt2Mat(path, "vecT");
	Matrix vecP = txt2Mat(path, "vecP");


	/*==========================================================================*/
	/*                 Part 1. Linear Least Square Regression                   */
	/*==========================================================================*/
	// Enter the temperature you want to find.
	double inputTemp1 = 150;

	Matrix coefficientZ = createMat(2, 1);
	coefficientZ = linearFit(vecT, vecP);
	
	double a0 = coefficientZ.at[0][0];
	double a1 = coefficientZ.at[1][0];
	double predict_P = (a0 * inputTemp1) + a1;
	
	printf("\n-------------------------------------------------------------");
	printf("\n         Part 1. Linear Least Square Regression              ");
	printf("\n-------------------------------------------------------------\n");
	printMat(coefficientZ, "Coefficient [a0] and a[1]");

	printf("Equation Form : f(x) = %fx + %f \n\n", a0, a1);
	printf("At Input Temperature %f, Pressure is %f \n\n\n", inputTemp1, predict_P);



	/*==========================================================================*/
	/*                  Part 2. Linear Spline Interpolation                     */
	/*==========================================================================*/
	double T_array[] = { 0, 5, 15, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100 };
	Matrix data_y = createMat(21, 1);
	Matrix T = createMat(21, 1);
	T = arr2Mat(T_array, 21, 1);

	// Enter the temperature you want to find.
	double inputTemp2 = 75;

	data_y = linearInterp(vecT, vecP, T);

	double find_P = 0;
	for (int i = 0; i < T.rows; i++)
	{
		if (inputTemp2 == T.at[i][0])
		{
			find_P = data_y.at[i][0];
		}
	}

	printf("\n-------------------------------------------------------------");
	printf("\n            Part 2. Linear Spline Interpolation              ");
	printf("\n-------------------------------------------------------------\n");
	printMat(data_y, "pressure :: y data");
	printf("At the input Temperature %f, Pressure is %f \n\n\n", inputTemp2, find_P);



	/*==========================================================================*/
    /*                       Deallocate memory							        */
    /*==========================================================================*/
    freeMat(vecT);			  freeMat(vecP);
	freeMat(coefficientZ);    freeMat(data_y);		freeMat(T);

	system("pause");
	return 0;
}