/*------------------------------------------------------------------------------------------
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Young-Keun Kim
Created         : 07-03-2019
Modified        : 03-01-2021
Language/ver    : C /  MSVS2017
Course          : Numerical method 2021-Spring
Description     : Tutorial of Bisection Method
/------------------------------------------------------------------------------------------*/

#include "myNM.h"
#include <math.h>

void main() {

    /*==========================================================================*/
    /*                     Solve by Newton-Raphson Method                       */
    /*==========================================================================*/

    /************      Variables declaration & initialization      ************/

    double _x = 1.4;               // try to choose 'x' close true value
    double _tol = 0.000001;
    double NRM_result;

    printf("------------------------------------------------------------------------------------\n");
    printf("         Newton-Raphson Method Results             \n");
    printf("------------------------------------------------------------------------------------\n");

    printf("Newton-Raphson Method:\n");
    NRM_result = newtonRaphson_2(_x, _tol);

    printf("Final Solution: %f \t", NRM_result);
    printf("\n \n \n");


    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    
    /*==========================================================================*/
    /*                  Assignment -     hybrid method                          */
    /*==========================================================================*/

    /************      Variables declaration & initialization      ************/
    
    double _x1 = 1.4;
    double Hybrid_result;

    printf("------------------------------------------------------------------------------------\n");
    printf("           Hybrid Method Results             \n");
    printf("------------------------------------------------------------------------------------\n");

    printf("Hybrid Method:\n");
    Hybrid_result = Hybrid(_x1, _tol);

    printf("Final Solution: %f \t", Hybrid_result);
    printf("\n \n \n");
}

