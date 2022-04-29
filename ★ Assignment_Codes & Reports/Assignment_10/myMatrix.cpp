/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Hong Se Hyun
Created          : 18-05-2021
Modified         : 25-06-2021
Language/ver     : C++ in MSVS2019

Description      : myMatrix.cpp
----------------------------------------------------------------*/

#include "myMatrix.h"

// Create Matrix with specified size
Matrix	createMat(int _rows, int _cols)
{
	// check matrix dimension
	if (_rows < 0 || _cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'createMat' function");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}		

	Matrix Out;
	// 1. Allocate row array first
	Out.at = (double**)malloc(sizeof(double*) * 2 * _rows);
	// 2. Then, allocate column 
	for (int i = 0; i < _rows; i++)
		Out.at[i] = (double*)malloc(sizeof(double) * 2 *_cols);
	// 3. Initialize row & column values of a matrix
	Out.rows = _rows;
	Out.cols = _cols;

	return Out;
}


// Free a memory allocated matrix
void	freeMat(Matrix _A)
{
	// 1. Free allocated column memory
	for (int i = 0; i < _A.rows; i++)
		free(_A.at[i]);
	// 2. Free allocated row memory
	free(_A.at);
}


// Create a matrix from a text file
Matrix	txt2Mat(std::string _filePath, std::string _fileName)
{
	std::ifstream file;
	std::string temp_string, objFile = _filePath + _fileName + ".txt";
	int temp_int = 0, nRows = 0;

	file.open(objFile);
	if (!file.is_open()) {
		printf("\n*********************************************");
		printf("\n  Could not access file: 'txt2Mat' function");
		printf("\n*********************************************\n");
		return createMat(0, 0);
	}
	while (getline(file, temp_string, '\t'))
		temp_int++;
	file.close();

	file.open(objFile);
	while (getline(file, temp_string, '\n'))
		nRows++;
	file.close();

	int nCols = (temp_int - 1) / nRows + 1;
	Matrix Out = createMat(nRows, nCols);

	file.open(objFile);
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nCols; j++) {
			file >> temp_string;
			Out.at[i][j] = stof(temp_string);
		}
	file.close();

	return Out;
}


// Print matrix
void	printMat(Matrix _A, const char* _name)
{
	printf("%s =\n", _name);
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++)
			printf("%15.6f\t", _A.at[i][j]);
		printf("\n");
	}
	printf("\n");
}


// initialization of Matrix elements
void	initMat(Matrix _A, double _val)
{
	// add your code here
	// add your code here

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			_A.at[i][j] = _val;

	// add your code here
	// add your code here
}


// Create matrix of all zeros
Matrix	zeros(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);
	// add your code here
	// add your code here

	initMat(Out, 0);

	// add your code here
	// add your code here

	return Out;
}


// Create matrix of all one
Matrix	ones(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);
	// add your code here
	// add your code here

	initMat(Out, 1);

	// add your code here
	// add your code here

	return Out;
}


// Create Identity Matrix
Matrix	eye(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++) {
		for (int j = 0; j < _cols; j++) {
			if (i == j) {
				Out.at[i][j] = 1;
			}
			else {
				Out.at[i][j] = 0;
			}
		}
	}
	return Out;
}


// Copy Matrix
Matrix	copyMat(Matrix _A)
{
	Matrix Out = createMat(_A.rows, _A.cols);

	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++) {
			Out.at[i][j] = _A.at[i][j];
		}
	}

	return Out;
}


//Make transpose Matrix _ This Function is only for Square Matrix !!
//Matrix	transpose(Matrix _A)
//{
//	Matrix Out = createMat(_A.rows, _A.cols);
//
//	for (int i = 0; i < _A.rows; i++) {
//		for (int j = 0; j < _A.cols; j++) {
//			Out.at[i][j] = _A.at[j][i];
//		}
//	}
//
//	return 	Out;
//}


// Make transpose Matrix
Matrix	transpose(Matrix _A)
{
	Matrix Out = createMat(_A.cols, _A.rows);

	for (int i = 0; i < _A.cols; i++) {
		for (int j = 0; j < _A.rows; j++) {
			Out.at[i][j] = _A.at[j][i];
		}
	}

	return 	Out;
}


// Norm
double	Norm(Matrix _c)
{
	double Out = 0;
	
	for (int i = 0; i < _c.rows; i++)
	{
		Out += (_c.at[i][0] * _c.at[i][0]);
	}
	Out = sqrt(Out);

	return Out;
}


// Multiply Matrix A*B and then save the result at C
void   multiply_Mat(Matrix _X, Matrix _Y, Matrix _Z)
{
	initMat(_Z, 0);
	
	for (int i = 0; i < _X.rows; i++) {
		for (int j = 0; j < _Y.cols; j++) {
			for (int k = 0; k < _X.cols; k++)
				_Z.at[i][j] += _X.at[i][k] * _Y.at[k][j];
		}
	}
}


// Find Max Value from Column Vector
double	max(Matrix _A)
{
	double max = _A.at[0][0];

	for (int i = 0; i < _A.rows ; i++)
	{
		if (max <= _A.at[i][0])
		{
			max = _A.at[i][0];
		}
	}
	
	return max;
}


// Find Min Value from Column Vector
double	min(Matrix _A)
{
	double min = _A.at[0][0];

	for (int i = 0; i < _A.rows; i++)
	{
		if (_A.at[i][0] != 0)
		{
			if (min >= _A.at[i][0])
			{
				min = _A.at[i][0];
			}
		}
	}

	return min;
}