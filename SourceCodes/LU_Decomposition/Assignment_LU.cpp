/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 22-03-2021
Language/ver     : C++ in MSVS2019

Description      : NM_MainTemplate.cpp
-------------------------------------------------------------------------------*/

#define Assignment	3		// enter your assignment number
#define eval		0		// set

#include "../../include/myNM.h"
//#include "myNM.h"

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
	/*   - You can change the variable names								4*/
	/*   - However, you must use the specified file name						*/
	/*	   : For each assignment, the file name will be notified on HISNET		*/
	/*==========================================================================*/

	Matrix matA = txt2Mat(path, "prob1_matA");
	Matrix vecb = txt2Mat(path, "prob1_vecb");


	// ERROR CHECKING 1
	if (matA.rows != matA.cols) {
		printf("\n*********************************************");
		printf("\n    ! ERROR ! It is not a square Matrix    ");
		printf("\n*********************************************\n");
		printf("\n \n \n \n \n");
	}

	// ERROR CHECKING 2
	if (vecb.cols != 1) {
		printf("\n*********************************************");
		printf("\n! ERROR ! input Matrix _b is not a column Vector ");
		printf("\n*********************************************\n");
		printf("\n \n \n \n \n");
	}


	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/

	Matrix vecd = createMat(vecb.rows, vecb.cols);
	Matrix matY = createMat(matA.rows, matA.cols);
	Matrix matP = createMat(matA.rows, matA.cols);
	Matrix matL = createMat(matA.rows, matA.cols);
	Matrix matU = createMat(matA.rows, matA.cols);
	Matrix x = createMat(matA.rows, vecb.cols);
	Matrix matX = createMat(matA.cols, vecb.cols);
	Matrix Ainv = createMat(matA.rows, matA.cols);

	Matrix identity = createMat(matA.rows, matA.cols);
	Matrix SP = createMat(matA.rows, 1);
	Matrix temp = createMat(1, matA.cols);
	Ainv = inv(matA);


	printf("\n-------------------------------------------------------------");
	printf("\n                      Input Matrix A                         ");
	printf("\n-------------------------------------------------------------\n");
	printMat(matA, "[ matA ]");
	

	/*==========================================================================*/
	/*				    	Print PLU Decomposition Result	      				*/
	/*==========================================================================*/
	LUdecomp(matA, matL, matU, matP);

	printf("\n-------------------------------------------------------------");
	printf("\n             P L U Decomposition of Matrix A                 ");
	printf("\n-------------------------------------------------------------\n");
	printMat(matP, "[ P ]");
	printMat(matL, "[ L ]");
	printMat(matU, "[ U ]");


	/*==========================================================================*/
	/*				  	Print the Solution by PLU Decomposition   				*/
	/*==========================================================================*/
	solveLU(matL, matU, matP, vecb, x);
	matX = matXvec(Ainv, vecb);
	printf("\n-------------------------------------------------------------");
	printf("\n            Solving matrix A by PLU Decomposition            ");
	printf("\n-------------------------------------------------------------\n");
	printMat(x, " X :: [ Using PLU Decomposition ] ");


	/*==========================================================================*/
	/*						Print the Solution by Inverse						*/
	/*==========================================================================*/
	printf("\n-------------------------------------------------------------");
	printf("\n   Inverse Matrix of A   &   Solving matrix A by Inverse     ");
	printf("\n-------------------------------------------------------------\n");
	printMat(Ainv, "Ainv");
	printMat(matX, " X :: [ Using Inverse Matrix Ainv ] ");


	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/

	freeMat(matA);		freeMat(vecb);		freeMat(vecd);		freeMat(matY);
	freeMat(matP);		freeMat(matL);		freeMat(matU);		freeMat(x);
	freeMat(matX);		freeMat(Ainv);		freeMat(identity);	freeMat(SP);		freeMat(temp);

	system("pause");
	return 0;
}