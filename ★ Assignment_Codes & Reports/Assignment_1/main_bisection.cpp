/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Young-Keun Kim
Created         : 07-03-2019
Modified        : 03-01-2021
Language/ver    : C /  MSVS2017 
Course          : Numerical method 2021-Spring
Description     : Tutorial of Bisection Method
/------------------------------------------------------------------------------------------*/

#include "myNM.h"

double bisectionNL(float _a0, float _b0, float _tol);  // Move this function to myNM.h

void main() {

    /*==========================================================================*/
    /*               Tutorial -     Bisection Method                            */
    /*==========================================================================*/

    /************      Variables declaration & initialization      ************/

    float tol = 0.00001;
    float a0 = 2.2; //need to change initial value
    float b0 = 3; //need to change initial value
    double BM_result;


    /************      Test NM Functions & Show Output            ************/
    printf("------------------------------------------------------------------------------------\n");
    printf("         Bisection Method Results             \n");
    printf("------------------------------------------------------------------------------------\n");

    printf("Bisection Method:\n");
    BM_result = bisectionNL(a0, b0, tol);

    printf("Final Solution: %f \t", BM_result);
    printf("\n");

    system("pause");
}


double func(double _x) // Move this function to myNM.c
{
    double F = 8 - 4.5*(_x - sin(_x));
    return F;
}



/* Bisection Method
   _a      : initial value #1
   _b      : initial value #2
   _tol   : tolerance   */

double bisectionNL(float _a0, float _b0, float _tol) // Move this function to myNM.c
{
    int k = 0;
    int Nmax = 1000;
    float a = _a0; 
    float b = _b0;
    float xn = 0;
    float ep = 1000;
    
    // You should check     if(func(a) * func(b) >= 0)
    // What should you do for this condition?

    do {
        xn = (a + b) / 2;
        ep = fabs(func(xn));
        printf("Iteration:%d \t", k);
        printf("X(n): %f \t", xn);
        printf("Tolerance: %.10f\n", ep);
         
        if (func(a) * func(xn) < 0)
            b = xn;
        else
            a = xn;
        k++;
    }while(k<Nmax && ep>_tol);

    
    return xn;
}

