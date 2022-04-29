/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Young-Keun Kim
Created         : 07-03-2019
Modified        : 03-01-2021
Language/ver    : C /  MSVS2017
Course          : Numerical method 2021-Spring
Description     : Tutorial of Bisection Method
/------------------------------------------------------------------------------------------*/

#include "../../include/myNM.h"
//#include "myNM.h"

void main() {

    /*==========================================================================*/
    /*                Assignment -     Bisection Method                         */
    /*==========================================================================*/

    /************      Variables declaration & initialization      ************/

    double _tol = 0.000001;
    double _a = 1.0;               // initial value #1
    double _b = 3.0;               // initial value #2
    double BM_result;

    printf("------------------------------------------------------------------------------------\n");
    printf("         Bisection Method Results             \n");
    printf("------------------------------------------------------------------------------------\n");

    printf("Bisection Method:\n");
    BM_result = bisectionNL(_a, _b, _tol);

    printf("Final Solution: %f \t", BM_result);
    printf("\n \n \n");


    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////


    /*==========================================================================*/
    /*                Assignment -     Newton-Raphson Method                    */
    /*==========================================================================*/

    /************      Variables declaration & initialization      ************/

    double _x = 1.0;               // try to choose 'x' close true value
    double NRM_result;

    printf("------------------------------------------------------------------------------------\n");
    printf("         Newton-Raphson Method Results             \n");
    printf("------------------------------------------------------------------------------------\n");
    
    printf("Newton-Raphson Method:\n");
    NRM_result = newtonRaphson(_x, _tol);

    printf("Final Solution: %f \t", NRM_result);
    printf("\n");

    system("pause");
}
