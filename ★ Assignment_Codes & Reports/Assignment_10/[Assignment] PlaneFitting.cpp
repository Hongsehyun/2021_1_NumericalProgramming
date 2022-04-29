/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Hong Se Hyun
Created          : 24-06-2021
Modified         : 25-06-2021
Language/ver     : C++ in MSVS2019

Description      : PlaneFit
----------------------------------------------------------------*/

#define Assignment	10		// enter your assignment number
#define PI		    3.14159265358979323846264338327950288

#include "myNM.h"
#include <math.h>
#include <stdio.h>


int main(int argc, char* argv[])
{
	/*	 [¡Ø DO NOT EDIT IT !!!]   Resources file path setting for evaluation	*/
	std::string path = "C:/NM_data_2021/Assignment" + std::to_string(Assignment) + "/";

#if eval
	path += "eval/";
#endif

	Matrix X_1_data = txt2Mat(path, "X_1_data");
	Matrix Y_1_data = txt2Mat(path, "Y_1_data");
	Matrix Z_1_data = txt2Mat(path, "Z_1_data");

	Matrix coefficientZ = createMat(3, 1);
	Matrix coefficientZ_y = createMat(3, 1);
	coefficientZ = planeFit(X_1_data, Y_1_data, Z_1_data);

	double a0 = coefficientZ.at[0][0];
	double a1 = coefficientZ.at[1][0];
	double a2 = coefficientZ.at[2][0];

	printf("\n-------------------------------------------------------------");
	printf("\n                      Part 1. Plane Fit                      ");
	printf("\n-------------------------------------------------------------\n");

	printMat(coefficientZ, "Coefficient a[0] a[1] a[2]");
	printf("Plane Equation : z = %f + %fy + %fx \n\n", a2, a1, a0);



	printf("\n--------------------------------------------------------------------------------------");
	printf("\n      Replace the plane equation for z with the plane equation for y as follows.      ");
	printf("\n--------------------------------------------------------------------------------------\n");
	coefficientZ_y.at[0][0] = a2 / (-a1);
	coefficientZ_y.at[1][0] = -1 / (-a1);
	coefficientZ_y.at[2][0] = a0 / (-a1);

	printf("Coefficients of Plane Equations for y\n\n");
	printMat(coefficientZ_y, "Coefficient a[0] a[1] a[2]");
	
	printf("Plane Equation : y = %f + %fz + %fx \n\n", coefficientZ_y.at[0][0], coefficientZ_y.at[1][0], coefficientZ_y.at[2][0]);
	printf("Plane Equation : %fx + y + %fz + %f = 0\n\n", -coefficientZ_y.at[2][0], -coefficientZ_y.at[1][0], -coefficientZ_y.at[0][0]);
	


	printf("\n-------------------------------------------------------------");
	printf("\n                          Calculation                        ");
	printf("\n-------------------------------------------------------------\n");
	// Normalized normal vector
	double a = -coefficientZ_y.at[2][0];
	double b = 1;
	double c = -coefficientZ_y.at[1][0];
	double d = -coefficientZ_y.at[0][0];
	
	Matrix NN_vector = createMat(3, 1);
	Matrix abc		 = createMat(1, 3);
	Matrix abc_T	 = createMat(3, 1);
	double abc_T_Norm= 0;

	//Matrix temp		 = createMat(3, 1);

	abc.at[0][0] = a;
	abc.at[0][1] = b;
	abc.at[0][2] = c;

	abc_T		= transpose(abc);
	abc_T_Norm  = Norm(abc_T);
	
	for (int i = 0; i < NN_vector.rows; i++){
		NN_vector.at[i][0] = abc_T.at[i][0] / abc_T_Norm;
	}

	printMat(NN_vector, "The Normalized Normal Vector");


	double alpha		= 0;
	double beta			= 0;
	double gamma		= 0;
	double aiming_angle = 0;
	double tilt_angle   = 0;
	
	alpha	= (180 / PI) * atan(c / b);
	beta	= (180 / PI) * atan(a / c);
	gamma	= (180 / PI) * atan(b / a);
	
	aiming_angle = 90 - gamma;
	tilt_angle = alpha;

	printf("alpha        = %f\n", alpha);
	printf("beta         = %f\n", beta);
	printf("gamma        = %f\n", gamma);
	printf("aiming angle = %f\n", aiming_angle);
	printf("tilt angle   = %f\n", tilt_angle);



	printf("\n\n\n");
	system("pause");
	return 0;
}