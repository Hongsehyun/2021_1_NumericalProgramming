#include "myNM.h"
#include <math.h>


// Tutorial BisectionNL function
double Tutorial_func(double _x)
{
    double F = 8 - 4.5 * (_x - sin(_x));
    return F;
}


// Assignment function f(x)
double func(double _x)
{
    double L = 4;
    double E = 70 * pow(10, 9);
    double I = 52.9 * pow(10, -6);
    double w_0 = 20000;
 
    double F = (w_0 / (360 * L * E * I)) * (7 * pow(L, 4) - 30 * pow(L, 2) * pow(_x, 2) + 15 * pow(_x, 4));
    // 원함수에서 한 번 미분한 식을 입력해줌. 즉, 원함수의 기울기를 나타내는 함수이며 단위는 Pa의 단위 N/M^2에 맞춰서 식을 전개함.

    return F;
}


// Assignment function f'(x)
double dfunc(double _x)
{
    double L = 4;
    double E = 70 * pow(10, 9);
    double I = 52.9 * pow(10, -6);
    double w_0 = 20000;

    double dF = (w_0 / (360 * L * E * I)) * (-60 * pow(L, 2) * _x + 60 * pow(_x, 3));

    return dF;
}


// Assignment Bonus Point function f(x)
double func_2(double _x)
{
    double F = (1 / _x) - 2;
    return F;
}


// Assignment Bonus Point function f'(x)
double dfunc_2(double _x)
{
    double dF = (-1 / pow(_x, 2));
    return dF;
}


// Bisection Code
double bisectionNL(double _a0, double _b0, double _tol)
{
    int k = 0;
    int Nmax = 1000;
    double _a = _a0;
    double _b = _b0;
    double xn = 0;
    double ep = 1000;

    // if f(a) * f(b) is above 0, bisection method cannot running.
    // Therefore, we should check if(func(a) * func(b) >= 0)
    if (func(_a) * func(_b) >= 0)
    {
        printf("새로운 초기값 a와 b를 입력하세요. \n");
        goto EXIT;     // Immediately close the function bisectionNL
    }
   
    do {
        xn = (_a + _b) / 2;
        ep = fabs(func(xn));
        printf("Iteration:%d \t", k);
        printf("X(n): %f \t", xn);
        printf("Tolerance: %.10f\n", ep);

        if (func(_a) * func(xn) < 0)
            _b = xn;
        else
            _a = xn;
        k++;
       
    } while (k<Nmax && ep>_tol);
 
EXIT:
    return xn;
}


// newtonRaphson Code
double newtonRaphson(double _x0, double _tol)
{
    //Assume that f(x) is continuous function.
    //Variable Initializing
    int k = 0;
    int Nmax = 1000;
    double _x = _x0;
    double ep = 1000;
    
    // if f'(_x) is 0, newton Raphson method cannot running.
    // Therefore, we should check if(f'(_x) == 0) 
    if (dfunc(_x) == 0)
    {
        printf("새로운 초기값 x를 입력하세요. \n");
        goto EXIT;     // Immediately close the function newtonRaphson
    }

    //Repeat until(k < Nmax&& ep > _tol)
    do {
        _x = _x - (func(_x) / dfunc(_x));
        ep = fabs(func(_x));
     
        printf("Iteration:%d \t", k);
        printf("X(n): %f \t", _x);
        printf("Tolerance: %.10f\n", ep);

        k++;

    } while (k<Nmax && ep>_tol);

EXIT:
    return _x;
}


// newtonRaphson Code for Hybrid Problem
double newtonRaphson_2(double _x0, double _tol)
{
    //Assume that f(x) is continuous function.
    //Variable Initializing
    int k = 0;
    int Nmax = 10;     // Setting the appropriate Nmax value to protect overflow
    double _x = _x0;
    double ep = 1000;

    // if f'(_x) is 0, newton Raphson method cannot running.
    // Therefore, we should check if(f'(_x) == 0) 
    if (func_2(_x) == 0)
    {
        printf("새로운 초기값 x를 입력하세요. \n");
        goto EXIT;     // Immediately close the function newtonRaphson
    }

    //Repeat until(k < Nmax&& ep > _tol)
    do {
        _x = _x - (func_2(_x) / dfunc_2(_x));
        ep = fabs(func_2(_x));

        printf("Iteration:%d \t", k);
        printf("X(n): %f \t", _x);
        printf("Tolerance: %.10f\n", ep);

        k++;

    } while (k<Nmax && ep>_tol);

EXIT:
    return _x;
}


// Final Hybrid Code
double Hybrid(double _x0, double _tol)
{
    //Assume that f(x) is continuous function.
    //Variable Initializing
    int k = 0;
    int Nmax = 1000;
    double _x1 = _x0;
    double _x2 = 0;
    double mid = 0;
    double ep = 1000;

    //Repeat until(k < Nmax&& ep > _tol)
    do {
        // since given function(func_2) seems to be diverging in Newton-Raphson Method, we should use bisection method for func_2
        // To use bisection method, we should initailize two variables
        ep = fabs(func_2(_x1));
        _x2 = _x1;
        _x1 = _x1 - (func_2(_x1) / dfunc_2(_x1));

        if (_x1 < 0.1 || _x1 > 2.0)
            // if the result value of newton-raphson method is out of boundary, [0.5 1.4]
            // we can assume that func_2 is diverging.
            // because, in the ideal newton-raphson method, the resulf must be close to solution [0.5] step by step. Like 1.3  1.1  0.9  0.7  0.55  0.51 and so on.
            // Therefore, if (_x < 0.5 || _x > 1.4), we should use bisection method and prevent diverging.
            // Change the boudary :: 0.1 / 2.0  due to the announcement related to hybrid.  [2021.03.22. 21:42]
        {
            int num = 0;
            while (_x1 <= 0)
            {
                mid = (_x1 + _x2) / 2;
                _x1 = mid;
                num++;
            }

            mid = (_x1 + _x2) / 2;
            if (func_2(_x2) * func_2(mid) < 0)
                _x1 = mid;
            else
                _x2 = mid;
        }
        // if _x is not in the condition [under 0.5 or over 1.4],
        // then we don't need to use this if algorithm.
        // by using newton raphson method, we can easily find the solution
        printf("Iteration:%d \t", k);
        printf("X(n): %f \t", _x1);
        printf("Tolerance: %.10f\n", ep);
        k++;
    } while (k<Nmax && ep>_tol);
    return _x1;
}


////failed-hybrid code
//double Hybrid(double _x0, double _tol)
//{
//    int k = 0;
//    int nmax = 1000;
//    double _x1 = _x0;
//    double _x2 = 0;
//    double mid = 0;
//    double ep = 1000;
//
//    do {
//         ep = fabs(func_2(_x1));
//
//        _x2 = _x1;
//        _x1 = _x1 - (func_2(_x1) / dfunc_2(_x1));
//
//        if (_x1 < 0.1 || _x1 > 2.0)
//        {          
//            mid = (_x1 + _x2) / 2;
//
//            if (func_2(_x1) <= 0)
//                _x1 = mid;
//            if (func_2(_x2) * func_2(mid) < 0)
//                _x1 = mid;
//            else
//                _x2 = mid;
//        }
//        
//        printf("iteration:%d \t", k);
//        printf("x(n): %f \t", _x1);
//        printf("tolerance: %.10f\n", ep);
//
//        k++;
//    } while (k<nmax && ep>_tol);
//
//    return _x1;
//}