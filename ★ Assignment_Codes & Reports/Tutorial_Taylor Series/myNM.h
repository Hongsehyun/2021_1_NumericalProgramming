/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Jan Park
Created          : 05-03-2021
Modified         : 05-03-2021
Language/ver     : C in MSVS2019

Description      : myNM.h
/----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H
#define		PI		3.14159265358979323846264338327950288419716939937510582

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Factorial function
extern double factorial(double _x);

// Taylor series approximation for sin(x) using pre-defined functions (input unit: [rad])
extern double sinTaylor(double _x);

// Taylor series approximation for sin(x) using pre-defined functions (input unit: [deg])
extern double sindTaylor(double _x);

// Taylor series approximation for sin(x) without using pre-defined functions (input unit: [rad])
extern double sinTaylor2(double _x);

// Function that reduced the computation cost of sinTaylor2 (input unit: [rad])
extern double sinTaylor3(double _x);

// TA's code for sin(x) approximation using Taylor series (input unit: [rad])
extern double sinTaylor4(double _x);

#endif