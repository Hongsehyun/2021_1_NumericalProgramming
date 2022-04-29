#pragma once
#define	_MY_NM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Tutorial Function
extern double Tutorial_func(double _x);

//Given Function
extern double func(double _x);
extern double dfunc(double _x);
extern double func_2(double _x);
extern double dfunc_2(double _x);

// BisectionNL function
extern double bisectionNL(double _a0, double _b0, double _tol);

//NewtonRhapson function
extern double newtonRaphson(double _x0, double _tol);
extern double newtonRaphson_2(double _x0, double _tol);

//Hybrid function
extern double Hybrid(double _x0, double _tol);
