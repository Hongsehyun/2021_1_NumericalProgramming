/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Jan Park
Created          : 10-05-2021
Modified         : 10-05-2021
Language/ver     : C++ in MSVS2019

Description      : [Tutorial]Differentiation_student.cpp
-------------------------------------------------------------------------------*/

#define Assignment   6		   // enter your assignment number
#define eval		 0         // set

#include "myNM.h"



// Return the dy/dx results for the input data. (truncation error: O(h^2))
Matrix	gradient(Matrix _x, Matrix _y);


Matrix	gradient_1D(Matrix _x, Matrix _y);



// Define a function that defines the target equation.
double myFunc(const double x);

// Return the dy/dx results for the target equation. (truncation error: O(h^2))
Matrix	gradientFunc(double func(const double x), Matrix xin);




int main(int argc, char* argv[])
{
	/*    [¡Ø DO NOT EDIT IT !!!]   Resources file path setting for evaluation   */
	std::string path = "C:/NM_data_2021/Assignment" + std::to_string(Assignment) + "/";

#if eval
	path += "eval/";
#endif

	/*==========================================================================*/
	Matrix t = txt2Mat(path, "Q1_vect");
	Matrix pos = txt2Mat(path, "Q1_vecx");
	/*==========================================================================*/

	// PART 1
	printf("\n**************************************************");
	printf("\n|                     PART 1.                    |");
	printf("\n**************************************************\n");

	Matrix vel = gradient(t, pos);
	Matrix acc = gradient(t, vel);

	printMat(t, "t");
	printMat(pos, "pos");
	printMat(vel, "vel");
	printMat(acc, "acc");





	// PART 2
	printf("\n**************************************************");
	printf("\n|                     PART 2.                    |");
	printf("\n**************************************************\n");

	Matrix xin = txt2Mat("", "Q2_vecxin");
	Matrix dydx = gradientFunc(myFunc, xin);

	printMat(xin, "xin");
	printMat(dydx, "dydx");

	system("pause");
	return 0;
}


// Return the dy/dx results for the input data. (truncation error: O(h^2))
// Move this function to myNM.h and myNM.cpp
Matrix	gradient(Matrix _x, Matrix _y) {

	int m = _x.rows;
	int n = _y.rows;

	// check if (m==n)


	Matrix df = createMat(m, 1);


	//Assuming constant h
	double h = _x.at[1][0] - _x.at[0][0];


	//Assuming n>2
	if (n > 2)
	{
		// 1. 3-point FW
		// df = [ -3f(x) + 4f(x+1) - f(x+2) ] / 2h
		df.at[0][0] = (-3 * _y.at[0][0] + 4 * _y.at[1][0] - _y.at[2][0]) / (2 * h);

		// 2. 2-point central diff
		// df = f(x+1)-f(x-1) / 2h
		for (int i = 1; i < n - 1; i++)
			df.at[i][0] = (_y.at[i + 1][0] - _y.at[i - 1][0]) / (2 * h);

		// 3. 3-point BW
		// df = [ f(x-2)-4f(x-1)+3f(x) ] / 2h
		df.at[n - 1][0] = (_y.at[n - 3][0] - 4 * _y.at[n - 2][0] - _y.at[n - 1][0]) / (2 * h);
	}
	
	else
	{
		df.at[0][0] = 0;   // use 2-point forward
		df.at[1][0] = 1;   // use 3-point forward
	}

	return df;
}


/*
// Return the dy/dx results for the input data. (truncation error: O(h^2))
// Move this function to myNM.h and myNM.cpp
Matrix	gradient_1D(Matrix _x, Matrix _y) {

	int n = _x.rows;
	Matrix df = createMat(n, 1);


	//Assuming constant h
	double h = _x.at[1][0] - _x.at[0][0];

	//Assuming n>2

	// 1. 3-point FW
	// df = [ -3f(x) + 4f(x+1) - f(x+2) ] / 2h
	df.at[0][0] = 

	// 2. 2-point central diff
	// df = f(x+1)-f(x-1) / 2h

	// 3. 3-point BW
	// df = [ f(x-2)-4f(x-1)+3f(x) ] / 2h


	//return df;
}
*/







// Define a function that defines the target equation.
double myFunc(const double x) {
	return  x * x * x;
}

// Return the dy/dx results for the target equation. (truncation error: O(h^2))
// Move this function to myNM.h and myNM.cpp
Matrix	gradientFunc(double func(const double x), Matrix xin) {
	
	int n = xin.rows;
	Matrix y = createMat(n, 1);
	
	// _y.at[0][0]


	// define y[0] to y[n-1]
	for (int i = 0; i < n; i++)
		y.at[i][0] = func(xin.at[i][0]);

	// Use Gradient() Numerical differentiation
	return df = gradient(xin, y);
	//return df;
	//Matrix df = gradient(xin, y);
	
	// gradient()
	Matrix df = createMat(n, 1);
	return createMat(0, 0);

} 