/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Hong Se Hyun
Created          : 09-05-2021
Modified         : 09-05-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H

#include "myMatrix.h"

// Assignment 5
extern Matrix	linearFit(Matrix _x, Matrix _y);

extern Matrix	linearInterp(Matrix _x, Matrix _y, Matrix _xq);

extern Matrix	arr2Mat(double* _1Darray, int _rows, int _cols);

#endif