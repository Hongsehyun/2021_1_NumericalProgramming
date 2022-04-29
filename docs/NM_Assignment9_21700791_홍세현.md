[TOC]

------





# Assignment ①

## Non-Linear Solver

## BisectionNL()

The bisection method is a bracketing method that uses an interval which includes the true solution x_T in [a,b], where a, b are the upper & lower bound for the solution.

Here, we have to  assume or check that the function f(x) is

​	▪ continuous

​	▪ it has a solution within the interval [a, b]



It successively reduces the interval range with iterations based on the following simple concept.

 	▪  If f(x) is continuous in [a, b] and the solution x_T is between [a, b], 

​		then the sign of f(a) and f(b) are opposite or f(a)f(b)<0.




#### Algorithm

1. Initially, observe the zero-crossing f(x)=0 and check if it has a root

![1-1  Bisection Explain](https://user-images.githubusercontent.com/84533279/122436399-c11ec380-cfd3-11eb-9613-290556ea4d5d.JPG)

2. Then, Choose an initial interval [a, b] such that f(a)f(b)<0

![1-2  Bisection Explain](https://user-images.githubusercontent.com/84533279/122436417-c419b400-cfd3-11eb-916d-52f78a2c1cef.JPG)

3. Estimate the next solution as ![1-0](https://user-images.githubusercontent.com/84533279/122436445-c976fe80-cfd3-11eb-8dfc-9dea2508c142.JPG)

4. Check whether x_T is closer to a or b

   If f(a)f(x_N) < 0, then x_T is within [a, Xn]

   If f(a)f(x_N) > 0, then x_T is within [Xn, b]

5. Update the next estimation x_N(k+1) with the new range [ a, x_N(k) ] or [ x_N(k),  b ]

 ![1-2 _](https://user-images.githubusercontent.com/84533279/122436475-cf6cdf80-cfd3-11eb-886c-9f3376fe8656.JPG)

6. Check the tolerance or rel. error. 

![1-3_](https://user-images.githubusercontent.com/84533279/122436564-e4497300-cfd3-11eb-851a-89cf8f043635.jPG)

7. Repeat from step 3) until the tolerance condition is satisfied

![1-3  Bisection Explain](https://user-images.githubusercontent.com/84533279/122436576-e6abcd00-cfd3-11eb-8fbe-f3e9f9cc5b17.JPG)

```c
double bisectionNL(double _a0, double _b0, double _tol);
```

#### **Parameters**

- **_a0** = Boundary Condition [Start Point]
- **_b0** = Boundary Condition [End Point]
- **_tol** = Tolerance

* **xn** = Return Value, the solution of Bisection Method

  

#### Example code

```C
/************      Variables declaration & initialization      ************/

    double _tol = 0.000001;
    double _a = 1.0;               // initial value #1
    double _b = 3.0;               // initial value #2
    double BM_result;

    printf("----------------------------------------------\n");
    printf("         Bisection Method Results             \n");
    printf("----------------------------------------------\n");

    printf("Bisection Method:\n");
    BM_result = bisectionNL(_a, _b, _tol);

    printf("Final Solution: %f \t", BM_result);
    printf("\n \n \n");
```







## Newton Raphson()

The Newton Raphson method solves for the root of an equation of the form f(x)=0.

In contrast to bisection method, this uses the slope of the equation to converge to the solution at a much faster rate. This method requires the f(x) is

​	▪ continuous

​	▪ differentiable

​	▪ start at the initial point near the true solution

![1-4  Newton Explain](https://user-images.githubusercontent.com/84533279/122436620-ef040800-cfd3-11eb-94bb-92be0b195829.JPG)

```c
double newtonRaphson(double _x0, double _tol);
```

#### **Parameters**

- **_x0** = Initial Point [need to closed to root]
- **_tol** = Tolerance

* **xn** = Return Value, the solution of Bisection Method



#### Example code

```C
/************      Variables declaration & initialization      ************/

    double _x = 1.0;               // try to choose 'x' close true value
    double NRM_result;

    printf("---------------------------------------------------\n");
    printf("         Newton-Raphson Method Results             \n");
    printf("---------------------------------------------------\n");
    
    printf("Newton-Raphson Method:\n");
    NRM_result = newtonRaphson(_x, _tol);

    printf("Final Solution: %f \t", NRM_result);
    printf("\n");

```







## Hybrid()

Hybrid method is combined version of Bisection method and Newton-Raphson method. 



If Newton-Raphson gives the solution out of bounds, use bisection method for the next estimation of x(n+1).

If Newton-Raphson is not decreasing fast enough, or seems to be diverging, use bisection method for x(n+1).

Otherwise use Newton-Raphson

```c
double Hybrid(double _x0, double _tol);
```

#### **Parameters**

- **_x0** = Initial Point [need to closed to root]
- **_tol** = Tolerance

* **xn** = Return Value



```c
/************      Variables declaration & initialization      ************/

double _x1 = 1.4;
double Hybrid_result;

printf("---------------------------------------------\n");
printf("           Hybrid Method Results             \n");
printf("---------------------------------------------\n");

printf("Hybrid Method:\n");
Hybrid_result = Hybrid(_x1, _tol);

printf("Final Solution: %f \t", Hybrid_result);
printf("\n \n \n");
```







## func()

You can enter the function you want to solve.

```c
double func(double _x);
```

#### **Parameters**

- **_x0** = Variables







## dfunc()

You can enter the differential form of the function you want to solve.

```c
double dfunc(double _x);
```

#### **Parameters**

- **_x0** = Variables







------







# Assignment ②

# Linear Solver

## gaussElim()

Using  Gaussian Elimination, We can to solve a linear equation.

This function does not implement a partial pivoting function.



There are some conditions,

1. Input Matrix _A must be a square Matrix
2. input Vector Matrix _b must be a column Vector
3. input Matrix _A's Diagonal Components must not be zero

```c
void	gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d);
```

#### **Parameters**

- **_A** = Input Matrix A

- **_b** = Input Vector b

- **_U** = Result of Gaussian Elimination to matrix A

- **_d** = Result of Gaussian Elimination to vector b

  

#### Example code

```C

	/*==========================================================================*/
	/*					Variables declaration & initialization					*/
	/*--------------------------------------------------------------------------*/
	Matrix matA = txt2Mat(path, "prob1_matA");
	Matrix vecb = txt2Mat(path, "prob1_vecb");

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
```







------







# Assignment ③

## LUdecomp()

By using LU Decomposition, we can find the solution of Ax=b without finding Inverse Matrix.



#### Algorithm

![2-3  LU Decomp](https://user-images.githubusercontent.com/84533279/122436671-f9260680-cfd3-11eb-9c1d-12e37397fd67.jpg)



Scaled Partial Pivoting method was used to do Partial Pivoting.

![2-5 Partial Pivoting](https://user-images.githubusercontent.com/84533279/122436691-fe835100-cfd3-11eb-9aec-0d8294eb8565.jpg)

```c
void	LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P)
```

#### **Parameters**

- **_A** = Input Matrix

- **_L** = Lower Triangular Matrix

- **_U** = Upper Triangular Matrix

- **_P** = Permutation Matrix

  

#### Example code

```c
Matrix matP = createMat(matA.rows, matA.cols);
Matrix matL = createMat(matA.rows, matA.cols);
Matrix matU = createMat(matA.rows, matA.cols);

Matrix identity = createMat(matA.rows, matA.cols);
Matrix SP = createMat(matA.rows, 1);
Matrix temp = createMat(1, matA.cols);

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
```







## SolveLU()

After LU Decomposition, we can find the solution of Linear Equation.



#### Algorithm

![2-2 solve LU](https://user-images.githubusercontent.com/84533279/122436723-03e09b80-cfd4-11eb-8418-b7682ded32a6.JPG)

```c
void	solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x);
```

#### **Parameters**

- **L** = Lower Triangular Matrix

- **U** = Upper Triangular Matrix 

- **P** = Permutation Matrix

- **b** = Input Vector b

- **x** = Solution of Linear Equation

  

#### Example code

```c
/*==========================================================================*/
/*				  	Print the Solution by PLU Decomposition   				*/
/*==========================================================================*/
solveLU(matL, matU, matP, vecb, x);
matX = solveX(Ainv, vecb);
printf("\n-------------------------------------------------------------");
printf("\n            Solving matrix A by PLU Decomposition            ");
printf("\n-------------------------------------------------------------\n");
printMat(x, " X :: [ Using PLU Decomposition ] ");
```







## inv()

Finding the inverse Matrix of Input Matrix A

we can find inverse Matrix from LU Decomposition and Gauss Jordan Elimination.



#### Algorithm ①

![2-4 LU Decomp Sol](https://user-images.githubusercontent.com/84533279/122436742-080cb900-cfd4-11eb-8719-10b7b79b92ce.jpg)



#### Algorithm ②

![2-6 Gauss Jordan Sol](https://user-images.githubusercontent.com/84533279/122436753-0b07a980-cfd4-11eb-95c1-445669b8d28d.JPG)

In this code, used the Gauss Jordan method.

```c
Matrix	inv(Matrix _A);
```

#### **Parameters**

- **_A** = Input Matrix A

  

#### Example code

```c
/*==========================================================================*/
/*						Print the Solution by Inverse						*/
/*==========================================================================*/

printf("\n-------------------------------------------------------------");
printf("\n   Inverse Matrix of A   &   Solving matrix A by Inverse     ");
printf("\n-------------------------------------------------------------\n");
printMat(Ainv, "Ainv");
printMat(matX, " X :: [ Using Inverse Matrix Ainv ] ");
```







## SolveX()

After Inv(), we can find the solution of Linear Equation.

#### 

#### Algorithm

![2-7 SolveX](https://user-images.githubusercontent.com/84533279/122436771-0e029a00-cfd4-11eb-8bc4-dba99bef0ce8.JPG)

This Function implement Multiplying Matrix and Vector.

Thus this function is same as 'matXvec' Function.

```c
Matrix	solveX(Matrix _Ainv, Matrix vecb)
```

#### **Parameters**

- **_Ainv** = Lower Triangular Matrix

- **vecb** = Upper Triangular Matrix 

  

#### Example code

```c
/*==========================================================================*/
/*				  	Print the Solution by PLU Decomposition   				*/
/*==========================================================================*/
solveLU(matL, matU, matP, vecb, x);
matX = solveX(Ainv, vecb);
printf("\n-------------------------------------------------------------");
printf("\n            Solving matrix A by PLU Decomposition            ");
printf("\n-------------------------------------------------------------\n");
printMat(x, " X :: [ Using PLU Decomposition ] ");
```







## matXvec()

By using this function, we can find the solution of Linear Equation.

```c
Matrix	matXvec(Matrix _Ainv, Matrix vecb)
```

#### **Parameters**

- **_Ainv** = Lower Triangular Matrix
- **vecb** = Upper Triangular Matrix 







## gaussjordan()

Gauss elimination method was making the equivalent triangular matrix.

This can be modified to make the matrix as the reduced row-echelon form of diagonal matrix with ones in the pivot elements



#### Algorithm

 Step 1 - Make [A b] Matrix from Linear Equation Ax = b

 Step 2 - Implement Gaussian Elimination with Pivoting

 Step 3 - Implement Back Substitution

 Step 4 -  Make Diagonal terms to 1

 Then, the solution is on the right sides.

```c
Matrix	gaussJordan(Matrix _A, Matrix _b)
```

#### **Parameters**

- **_A** = Input Matrix A 
- **_b** = Input Vector b







------





# Assignment ④

## QRdecomp()

To get the eigenvalue of matrix A, we need to find Upper triangular Matrix and Orthogonal Matrix of A.

This process is called QR decomposition, and It can be written as follows.

![4-2 QR Decomp](https://user-images.githubusercontent.com/84533279/122436803-14911180-cfd4-11eb-9ca8-ba6f25b81ef5.JPG)



Then, how can we decompose A into Q, R?

   Method 1) Using Householder matrix

   Method 2) Gram-Schmidt factorization



In this function, we used Householder matrix

![5-3 Household Mat](https://user-images.githubusercontent.com/84533279/122436817-178c0200-cfd4-11eb-8d0e-304c10fc0cdf.JPG)

```c
void QRdecomp(Matrix _A, Matrix _H, Matrix _Q, Matrix _R)
```

#### **Parameters**

- **_A** = Input Matrix A
- **_H** = HouseHolder Matrix H
- **_Q** = Orthogonal Matrix Q
- **_R** = Upper Triangular Matrix R







## eig()

The eigenvalue of a square matrix A are defined as 

![5-5 eigenvalue](https://user-images.githubusercontent.com/84533279/122436838-1a86f280-cfd4-11eb-8e9d-3e404ade802a.JPG)



After QR Decomposition, we can easily find eigenvalue through the process below.

![5-2 QR Decomp](https://user-images.githubusercontent.com/84533279/122436842-1ce94c80-cfd4-11eb-9cc3-cea6dc10c4d4.JPG)

```c
Matrix	eig(Matrix _A, Matrix _H, Matrix _Q, Matrix _R)
```

#### **Parameters**

- **_A** = Input Matrix A
- **_H** = HouseHolder Matrix H
- **_Q** = Orthogonal Matrix Q
- **_R** = Upper Triangular Matrix R



#### Example code

```c
/*==========================================================================*/
/*                            Print EigenValue                          */
/*==========================================================================*/
// QR Decomposition Matrix Setting
    Matrix ATAmatH = createMat(ATA.rows, ATA.cols);
    Matrix ATAmatQ = createMat(ATA.rows, ATA.cols);
    Matrix ATAmatR = createMat(ATA.rows, ATA.cols);
    printf("\n--------------------------------------------------------------");
    printf("\n                EigenValues  of  A T A                        ");
    printf("\n--------------------------------------------------------------\n");
    eigenATA = eig(ATA, ATAmatH, ATAmatQ, ATAmatR);
    printMat(eigenATA, "eigenValue of ATA");
```







## condNumber()

the condition number of a function measures how much the output value of the function can change for a small change in the input argument.

This is used to measure how sensitive a function is to changes or errors in the input, and how much error in the output results from an error in the input. 



#### Algorithm

![5-4 Cond Num](https://user-images.githubusercontent.com/84533279/122436861-21156a00-cfd4-11eb-93f0-fe92b9f7c1c9.JPG)

```c
double	condNumber(Matrix _A)
```

#### **Parameters**

- **_A** =  Input Matrix A



#### Example code

```c
/*==========================================================================*/
/*                            Condition Number                              */
/*==========================================================================*/
printf("\n-------------------------------------------------------------");
printf("\n                 Condition Number of Mat A                   ");
printf("\n-------------------------------------------------------------\n");
Cond_NumberA = condNumber(matA);
printf("%f\t", Cond_NumberA);
printf("\n\n");
```







------







# Assignment ⑤

## linearFit()

Curve Fitting with a line (first degree polynomial) is finding the best parameters (a0, a1) for the linear equation

![6-1 Func](https://user-images.githubusercontent.com/84533279/122436872-2377c400-cfd4-11eb-8339-edf43dd18e8d.JPG)

that gives the least total error when a set of data (xi,yi) are given.

![6-2 Curve Fit](https://user-images.githubusercontent.com/84533279/122436884-2672b480-cfd4-11eb-86c7-3a76c79c1200.JPG)

Since the goal is to find a function with the least error, let's focus on the value of zero when differentiated,.

![6-3 Curve FIt](https://user-images.githubusercontent.com/84533279/122436893-28d50e80-cfd4-11eb-9751-2b5f94644a60.JPG)

```c
Matrix	linearFit(Matrix _x, Matrix _y)
```

#### **Parameters**

- **_x** = Input data set _x

- **_y** = Input data set _y

  

#### Example code

```C
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
```







## linearInterp()

In this Function, we used Lagrange Polynomial.

![6-4 Lagrange](https://user-images.githubusercontent.com/84533279/122436909-2bcfff00-cfd4-11eb-86dc-fbd3ec99b2e9.JPG)

note that, 10 data points returns 9 equations.

```
Matrix	linearInterp(Matrix _x, Matrix _y, Matrix _xq)
```

#### **Parameters**

- **_x** = Input data set _x

- **_y** = Input data set _y

- **xq** = query Input.

  

#### Example code

```c
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
```







------







# Assignment ⑥

## gradient()

This chapter will cover how to calculate a numerical approximation of the derivative from the given discrete datasets.

From this table, we can easily differentiate any datas or functions.

![7-1 Diff](https://user-images.githubusercontent.com/84533279/122436919-2ecaef80-cfd4-11eb-8152-6ca8e223509f.JPG)

```c
Matrix	gradient(Matrix _x, Matrix _y)
```

#### **Parameters**

- **_x** = Input variable x [ x in the form of dxdy ]

- **_y** = Input variable y [ y in the form of dxdy ]

  

#### Example code

```c
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
```







## gradient1D()

Differentiate use a one-dimensional array.

```C
void	gradient1D(double x[], double y[], double dydx[], int m) {
```

#### **Parameters**

- **x[]** = Input variable x [ x in the form of dxdy ]

- **y[]** = Input variable y [ y in the form of dxdy ]

- **dydx[]** = Store Results

- **m** = number of data

  

#### Example code

```C
// PART 1-2 :: gradient with Array
printf("\n**************************************************");
printf("\n|                     PART 1.2                   |");
printf("\n|                                                |");
printf("\n|             - Gradient with Array -            |");
rintf("\n**************************************************\n");
double t_arr[21] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6 ,3.8 ,4.0 };
double pos_arr[21] = { -5.87, -4.23, -2.55, -0.89, 0.67, 2.09, 3.31, 4.31, 5.06, 5.55, 5.78, 5.77, 5.52, 5.08, 4.46, 3.72, 2.88, 2.0, 1.1, 0.23, -0.59 };
double vel_arr[21] = { 0.0 };
double acc_arr[21] = { 0.0 };
	
//int m = sizeof(x_arr) / sizeof(double);
//int n = sizeof(y_arr) / sizeof(double);
	
int m = 21;
gradient1D(t_arr, pos_arr, vel_arr, m);
gradient1D(t_arr, vel_arr, acc_arr, m);

printf("vel_arr = \n");
printArr(vel_arr);
printf("acc_arr = \n");
printArr(acc_arr);
```







## gradientFunc()

Differentiate from a function.

```C
Matrix	gradientFunc(double myFunc(const double x), Matrix xin) {
```

#### **Parameters**

- **double myFunc(const double x)** = Input function

- **xin** = Input data set

  

#### Example code

```c
	// PART 2-1 :: gradient [ From a Function ]
	printf("\n**************************************************");
	printf("\n|                     PART 2.1                   |");
	printf("\n|                                                |");
	printf("\n|         - Gradient from the Function -         |");
	printf("\n**************************************************\n");

	Matrix dydx = gradientFunc(myFunc, xin);

	printMat(xin, "xin");
	printMat(dydx, "dydx");
```







## myFunc()

Define function

```c
double myFunc(const double x)
```

#### **Parameters**

- **_x** = Variable

  

#### Example code

```c
double myFunc(const double x) {
	double y = x * x * x;
	return y;
}
```







------







# Assignment ⑦

## IntegrateRect()

Rectangular Method for Discrete Data Input

Let's use the result of the integral is the area of the function.

![8-1 Rect](https://user-images.githubusercontent.com/84533279/122436944-34c0d080-cfd4-11eb-8371-b0ceae1db7e3.JPG)

```
double IntegrateRect(double _x[], double _y[], int _m) 
```

#### **Parameters**

- **x[]** = Input Dataset X

- **y[]** = Input Dataset X

- **_m** = Number of data

  

#### Example code

```c
double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
int M = sizeof(x) / sizeof(x[0]);

// Rectangular Method for Discrete Data Input
double I_rect = IntegrateRect(x, y, M);
printf("I_rect  = %f\n\n", I_rect);
```







## IntegralMid()

Mid Point Method for Discrete Data Input

Let's use the result of the integral is the area of the function.

![8-2 mid](https://user-images.githubusercontent.com/84533279/122436964-38545780-cfd4-11eb-8ec3-c0a4ae82423f.JPG)

```c
double IntegralMid(double _x[], double _y[], int _m)
```

#### **Parameters**

- **x[]** = Input Dataset X

- **y[]** = Input Dataset X

- **_m** = Number of data

  

#### Example code

```c
double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
int M = sizeof(x) / sizeof(x[0]);

// Mid Point Method for Discrete Data Input
double I_midpoint = IntegralMid(x, y, M);
printf("I_midpoint = %f\n\n", I_midpoint)
```







## trapz()

Trapezoidal Method for Discrete Data Input

Let's use the result of the integral is the area of the function.

![8-3 trap](https://user-images.githubusercontent.com/84533279/122436999-4013fc00-cfd4-11eb-8814-284bf1a38dd8.JPG)

![8-3 trap2](https://user-images.githubusercontent.com/84533279/122437002-40ac9280-cfd4-11eb-9ee7-5307a22d4887.JPG)


```c
double	trapz(double _x[], double _y[], int _m)
```

#### **Parameters**

- **x[]** = Input Dataset X

- **y[]** = Input Dataset X

- **_m** = Number of data

  

#### Example code

```c
double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
int M = sizeof(x) / sizeof(x[0]);

// Trapezoidal Method for Discrete Data Input
double I_trapz = trapz(x, y, M);
printf("I_trapz = %f\n\n", I_trapz);
```







## integral()

Simpson 1/3 Method for Function Input.

We need 3 points for unique quadratic polynomial

![8-4 simpson1](https://user-images.githubusercontent.com/84533279/122437017-44d8b000-cfd4-11eb-9a71-a539ed0f2277.JPG)

conditions : For N even numbers, same intervals

```c
double   integral(double func(const double _x), double _a, double _b, int _n) {
```

#### **Parameters**

- **double func(const double _x)** = Input function

- **_a** = Integration Start Point

- **_b** = Integration End Point

- **_n** = number of intervals

  

#### Example code

```c
double a = -1;		// 적분구간 a
double b = 1;		// 적분구간 b
double n = 12;      // Number of Intervals = 12

// Simpson 1/3 Method for Function Input
double I_simpson13 = integral(myFunc, a, b, n);
printf("I_simpson1/3  = %f\n\n", I_simpson13);
```







## integral38()

Simpson 3/8 Method for Function Input.

We need 4 points for unique cubic polynomial.

![8-5 simpson2](https://user-images.githubusercontent.com/84533279/122437032-4904cd80-cfd4-11eb-87e8-1013a8712bc6.JPG)

conditions : For N multiple of 3 and, same intervals

```c
double   integral38(double func(const double _x), double _a, double _b, int _n)
```

#### **Parameters**

- **double func(const double _x)** = Input function

- **_a** = Integration Start Point

- **_b** = Integration End Point

- **_n** = number of intervals

  

#### Example code

```c
double a = -1;		// 적분구간 a
double b = 1;		// 적분구간 b
double n = 12;      // Number of Intervals = 12

// Simpson 3/8 Method for Function Input
double I_simpson38 = integral38(myFunc, a, b, n);
printf("I_simpson3/8  = %f\n\n", I_simpson38);
```







------







# Assignment ⑧

## odeEU()

Euler Method

![9-1 Euler](https://user-images.githubusercontent.com/84533279/122437040-4b672780-cfd4-11eb-912a-09951ec84b61.JPG)

```C
void  odeEU(double func(const double t, const double v), double y[], double t0, double tf, double h, double y0)
```

#### **Parameters**

- **double func(const double t, const double v)** = slope function [returns dy/dt]

- **y[]** = 1-D array for y(t)

- **t0** = start point

- **tf** = finish point

- **h** = interval

- **y0** = Initial State

  

#### Example code

```C
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
```







## odeEM()

Euler Modified Method

![9-2 Euler Mod](https://user-images.githubusercontent.com/84533279/122437055-4efaae80-cfd4-11eb-9dff-d63fd8132dac.JPG)

```C
void odeEU_new(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0) 
```

#### **Parameters**

- **double func(const double t, const double v)** = slope function [returns dy/dt]

- **y[]** = 1-D array for y(t)

- **t0** = start point

- **tf** = finish point

- **h** = interval

- **y0** = Initial State

  

#### Example code

```C
// PART 2. Euler Modified Method
printf("\n\n\n\n**************************************************");
printf("\n           PART 2. Euler Modified Method          ");
printf("\n**************************************************\n");
// Initial Value
double value_EM[101] = { 0, };
/ Euler Modified Method
printf("Euler Modified Method :: Data of v\n\n");
odeEM(func, value_EM, t0, tf, interval);
```







------







# Assignment ⑨

## odeRK2()

It is a family of single-step, explicit NM solution for 1st ODE. The slope is obtained from a weighted sum of gradient at multiple points.

![10-1 RK](https://user-images.githubusercontent.com/84533279/122437072-53bf6280-cfd4-11eb-974c-97ac605e3ee2.JPG)

General Form of Runge - kutta 2 method uses two slopes to find next value of y.

And also, General Form of Runge - kutta 2 method set C1 = C2 = 0.5 and alpha = beta = 1. Then It is equal to Euler Modified method.

![10-2 RK](https://user-images.githubusercontent.com/84533279/122437084-56ba5300-cfd4-11eb-8e0c-f236193cc110.JPG)

```C
void odeRK2(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
```

#### **Parameters**

- **double odefunc(const double t, const double v)** = slope function [returns dy/dt]

- **y[]** = 1-D array for y(t)

- **t0** = start point

- **tf** = finish point

- **h** = interval

- **y0** = Initial State

  

#### Example code

```C
	//Parameter Definitions
	double a = 0;
	double b = 0.1;
	double h = 0.001;
	unsigned int N = (b - a) / h + 1;
	double y_EU[200] = { 0 };				//Cannot use y_EU[N]
	double y_EM[200] = { 0 };
	double y_RK2[200] = { 0 };
	double y_RK4[200] = { 0 };

	// Initial value
	double v0 = 0;

	// Exercise 1: Create a general form for RK2
	odeRK2(odeFunc_rc, y_RK2, a, b, h, v0);

	// Print outputs
	printf("/*-----------------------*/\n");
	printf("/ Single of 1st Order ODE /\n");
	printf("/*-----------------------*/\n");
	printf(" - Total number of data N=%d \n", N);
	for (int i = 0; i < N; i++)
		printf("%f\tyEU= %f\tyEM= %f\tyRK2= %f\tyRK4= %f\n", a + i * h, y_EU[i], y_EM[i], y_RK2[i], y_RK4[i]);
	printf("\n");
```







## odeRK4()

General Form of Runge - kutta 4 method uses 4 slopes to find next value of y.

And also General Form of Runge - kutta 4 method set C1 = C4 = 1/6 , C2 = C3 = 2 / 6 and alpha = beta = 1/2.

![10-3 RK](https://user-images.githubusercontent.com/84533279/122437092-5a4dda00-cfd4-11eb-9fec-45ab49895025.jpg)

```C
void odeRK4(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
```

#### **Parameters**

- **double odefunc(const double t, const double v)** = slope function [returns dy/dt]

- **y[]** = 1-D array for y(t)

- **t0** = start point

- **tf** = finish point

- *h** = interval

- **y0** = Initial State

  

#### Example code

```C
	//Parameter Definitions
	double a = 0;
	double b = 0.1;
	double h = 0.001;
	unsigned int N = (b - a) / h + 1;
	double y_EU[200] = { 0 };				//Cannot use y_EU[N]
	double y_EM[200] = { 0 };
	double y_RK2[200] = { 0 };
	double y_RK4[200] = { 0 };

	// Initial value
	double v0 = 0;

	// Exercise 2: Create the standard form  for RK4
	odeRK4(odeFunc_rc, y_RK4, a, b, h, v0);

	// Print outputs
	printf("/*-----------------------*/\n");
	printf("/ Single of 1st Order ODE /\n");
	printf("/*-----------------------*/\n");
	printf(" - Total number of data N=%d \n", N);
	for (int i = 0; i < N; i++)
		printf("%f\tyEU= %f\tyEM= %f\tyRK2= %f\tyRK4= %f\n", a + i * h, y_EU[i], y_EM[i], y_RK2[i], y_RK4[i]);
	printf("\n");
```







## sys2RK2()

Applying Runke-Kutta method to 2nd order ODE.

A second order ODE with dependent variable y(t) and one independent variable t.

Here, I introduce one tip for approaching quadratic differential equations.

This Equation is from m-c-k system.

![10-4 2nd order RK explain](https://user-images.githubusercontent.com/84533279/122437110-5d48ca80-cfd4-11eb-9432-40a91cc0ff20.JPG)

```C
void sys2RK2(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init)
```

#### **Parameters**

- **void odeFunc_sys2(const double t, const double Y[], double dYdt[])** = slope function [returns dy/dt and dz/dt]

- **y1[]** = 1-D array for y(t)

- **y2[]** = 1-D array for z(t) [ z(t) = dydt(t) ]

- **t0** = start point

- **tf** = finish point

- **h** = Interval

- **y1_init** = Initial State

- **y2_init** = Initial State

  

#### Example code

```c
//Parameter Definitions
double t0 = 0;
double tf = 1;
h = 0.01;
N = (tf - t0) / h + 1;
		
double y_sysRK2[200] = { 0 };
double v_sysRK2[200] = { 0 };
	
// Initial values
double y0 = 0;
v0 = 0.2;

// ODE solver: RK2
sys2RK2(odeFunc_mck, y_sysRK2, v_sysRK2, t0, tf, h, y0, v0);
```







## sys2RK4()

Applying Runke-Kutta method to 2nd order ODE.

```c
void sys2RK4(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init) {
```

#### **Parameters**

- **void odeFunc_sys2(const double t, const double Y[], double dYdt[])** = slope function [returns dy/dt and dz/dt]

- **y1[]** = 1-D array for y(t)

- **y2[]** = 1-D array for z(t) [ z(t) = dydt(t) ]

- **t0** = start point

- **tf** = finish point

- **h** = Interval

- **y1_init** = Initial State

- **y2_init** = Initial State

  

#### Example code

```c
//Parameter Definitions
double t0 = 0;
double tf = 1;
h = 0.01;
N = (tf - t0) / h + 1;
		
double y_sysRK4[200] = { 0 };
double v_sysRK4[200] = { 0 };

// Initial values
double y0 = 0;
v0 = 0.2;

// Create the standard form  for RK4 for 2nd order	
sys2RK4(odeFunc_mck, y_sysRK4, v_sysRK4, t0, tf, h, y0, v0);
```







------







# Matrix Functions





## addMat()

Matrix size of Input Matrix A and B must be equal.

This function performs mat A + mat B.

```c
Matrix	addMat(Matrix _A, Matrix _B)
```

#### **Parameters**

- **_A** = Input Matrix A
- **_B** = Input Matrix B







## fwdSub()

A function that changes the form of given matrix  such as:

![2-8 forward sub Form](https://user-images.githubusercontent.com/84533279/122437147-65086f00-cfd4-11eb-8e54-157ea7cb0cf8.JPG)

```c
void	fwdSub(Matrix _L, Matrix _d)
```

#### **Parameters**

- **_L** = Lower Triangular Matrix [ Result of forward substitution to matrix _A ]
- **_d** = Result of forward substitution to Vector Matrix _b







## backSub()

A function that changes the form of given matrix  such as:

![2-9 Back sub form](https://user-images.githubusercontent.com/84533279/122437154-676ac900-cfd4-11eb-9e85-ca4c2b23d648.JPG)

#### Algorithm

![2-1  BackSub](https://user-images.githubusercontent.com/84533279/122437162-69cd2300-cfd4-11eb-9a05-97b6bc90b134.JPG)

```c
void	backSub(Matrix _U, Matrix _d, Matrix _x)
```

#### **Parameters**

- **_U** = Upper Triangular Matrix [ Result of Gaussian Elimination to matrix _A ]
- **_d** = Result of Gaussian Elimination to Vector Matrix _b

* **x** = Solution of Linear Equation







## initMat()

 initialization of Matrix elements

```c
void	initMat(Matrix _A, double _val)
```

#### **Parameters**

- **_A** = Input Matrix
- **_val** = Reset Value







## zeros()

Create a matrix with all values equal to zero

```c
Matrix	zeros(int _rows, int _cols)
```

#### **Parameters**

- **_rows** = Number of rows of matrices you want to create
- **_cols** = Number of cols of matrices you want to create







## ones()

Create a matrix with all values equal to one

```c
Matrix	ones(int _rows, int _cols)
```

#### Parameters

- **_rows** = Number of rows of matrices you want to create
- **_cols** = Number of cols of matrices you want to create







## eye()

Create a identity Matrix.

Identity matrix is a matrix in which only the main diagonal component has a value of 1 and all the remaining components have a value of 0.

![2-10 Identity Mat](https://user-images.githubusercontent.com/84533279/122437186-6e91d700-cfd4-11eb-8266-eb967292aa77.JPG)

```c
Matrix	eye(int _rows, int _cols)
```

#### Parameters

- **_rows** = Number of rows of matrices you want to create
- **_cols** = Number of cols of matrices you want to create







## copyMat()

Copy Matrix.

```c
Matrix	copyMat(Matrix _A)
```

#### Parameters

- **_A** = Input Matrix







## transpose()

Function that returns the transposition matrix.

The transpose of a matrix is an operator which flips a matrix over its diagonal

```c
Matrix	transpose(Matrix _A)
```

#### Parameters

- **_A** = Input Matrix







## Norm()

Function that calculates the Norm

![4-1 Norm](https://user-images.githubusercontent.com/84533279/122437205-73568b00-cfd4-11eb-8cb2-15b1560ffe0c.JPG)

In this function, we used 2-norm to calculate the condition number:

```c
double	Norm(Matrix _c)
```

#### Parameters

- **_c** = Input Matrix c







## multiply_Mat()

Multiply Matrix and then save the result at _Z

In this function, mat _Z serves as a temporary repository.

```c
void   multiply_Mat(Matrix _X, Matrix _Y, Matrix _Z)
```

#### Parameters

- **_X** = Input Matrix X
- **_Y** = Input Matrix Y
- **_Z** = The result of multiplication of X and Y







## max()

Find Max Value from Column Vector

```c
double	max(Matrix _A)
```

#### Parameters

- **_A** = Input Vector







## min()

Find Min Value from Column Vector.

```c
double	min(Matrix _A)
```

#### Parameters

- **_A** = Input Vector
