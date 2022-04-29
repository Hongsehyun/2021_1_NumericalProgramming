/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 22-03-2021
Language/ver     : C++ in MSVS2019

Description      : NM_MainTemplate.cpp
-------------------------------------------------------------------------------*/

#define Assignment	3		// enter your assignment number
#define eval		0		// set 0

#include "myNM.h"

int main(int argc, char* argv[])
{
	/*	 [¡Ø DO NOT EDIT IT !!!]   Resources file path setting for evaluation	*/
	std::string path = "C:/NM_data_2021/Assignment" + std::to_string(Assignment) + "/";

#if eval
	path += "eval/";
#endif

	/*==========================================================================*/
	/*					Variables declaration & initialization					*/
	/*--------------------------------------------------------------------------*/
	/*   - You can change the variable names									*/
	/*   - However, you must use the specified file name						*/
	/*	   : For each assignment, the file name will be notified on HISNET		*/
	/*==========================================================================*/

	Matrix matA = txt2Mat(path, "prob1_matA");
	Matrix vecb = txt2Mat(path, "prob1_vecb");


	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/

	Matrix matU = createMat(matA.rows, matA.cols);
	Matrix vecd = createMat(vecb.rows, vecb.cols);
	Matrix matX = createMat(matA.cols, vecb.cols);

	printMat(matA, "matA");
	printMat(vecb, "vecb");

	gaussElim(matA, vecb, matU, vecd);
	backSub(matU, vecd, matX);
	
	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/

	printMat(matU, "matU");
	printMat(vecd, "vecd");
	printMat(matX, "matX");

	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	freeMat(matA);		freeMat(vecb);
	freeMat(matU);		freeMat(vecd);

	system("pause");
	return 0;
}