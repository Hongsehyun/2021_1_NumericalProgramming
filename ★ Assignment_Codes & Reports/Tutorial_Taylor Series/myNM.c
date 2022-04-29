/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Jan Park
Created          : 05-03-2021
Modified         : 05-03-2021
Language/ver     : C in MSVS2019

Description      : myNM.c
/----------------------------------------------------------------*/

#include "myNM.h"

// factorial function
double factorial(double _x)
{
	if (_x <= 1)
		return 1;
	else
		return _x * factorial(_x - 1);
}

//  Taylor series approximation for sin(x) using pre-defined functions (input unit: [rad])
double sinTaylor(double _x) 
{
	int N_max = 20;
	double epsilon = 1e-5;

	double S_N = 0, S_N_prev = 0, rel_chg = 0;
	int N = 0;

	do {
		N++;
		S_N_prev = S_N;
		S_N = 0;
		for (int k = 0; k < N; k++)
			S_N += pow(-1, k) * pow(_x, 2 * k + 1) / factorial(2 * k + 1);

		rel_chg = fabs((S_N - S_N_prev) / S_N_prev);

	} while (N < N_max && rel_chg >= epsilon);

	return S_N;
}

// Taylor series approximation for sin(x) using pre-defined functions (input unit: [deg])
double sindTaylor(double _x) 
{
	return sinTaylor(_x * PI / 180);
}

// Taylor series approximation for sin(x) without using pre-defined functions (input unit: [rad])
double sinTaylor2(double _x) 
{
	int N_max = 20;
	double epsilon = 1e-5;

	double S_N = 0, S_N_prev = 0, rel_chg = 0;
	int N = 0;

	do {
		N++;
		S_N_prev = S_N;
		S_N = 0;

		for (int k = 0; k < N; k++)
		{
			int sign_part = 1;
			for (int i = 1; i <= k; i++)
				sign_part *= -1;

			double pow_part = 1;
			for (int i = 1; i <= 2 * k + 1; i++)
				pow_part *= _x;

			double fac_part = 1;
			for (int i = 1; i <= 2 * k + 1; i++)
				fac_part *= i;

			S_N += sign_part * pow_part / fac_part;
		}

		rel_chg = fabs((S_N - S_N_prev) / S_N_prev);

	} while (N < N_max && rel_chg >= epsilon);

	return S_N;
}

// Function that reduced the computation cost of sinTaylor2 (input unit: [rad])
double sinTaylor3(double _x) 
{
	int N_max = 20;
	double epsilon = 1e-5;

	double S_N = 0, S_N_prev = 0, rel_chg = 0;
	int N = 0;
		   
	int sign_part = -1;
	double pow_part = 1 / _x;
	double fac_part = 1;
	
	do {
		N++;
		S_N_prev = S_N;

		sign_part *= -1;
		pow_part *= _x * _x;
		fac_part = max(fac_part * (2 * N - 2) * (2 * N - 1), 1); // max 함수는 입력된 인자 중에 가장 큰 값만 출력해줌.
																 // 여기서 max를 쓴 이유는, 처음 실행되는 단계인 N=1인 상태에서 fac_part가 0이 나오기 때문에,
																 // 이 경우를 보완하고자 max(0,1)값을 통해 1이 강제적으로 출력되게 만든 것!
		S_N += (sign_part * pow_part / fac_part);
		
		rel_chg = fabs((S_N - S_N_prev) / S_N_prev);
		
	} while (N < N_max && rel_chg >= epsilon);

	return S_N;
}

// Function to prevent overflow of sinTaylor3 (input unit: [rad]) - 이 부분 역시 SinTaylor3의 코드를 단순화, 간략화하여 연산 속도를 늘려줌
double sinTaylor4(double _x) 
{
	int N_max = 20;
	double epsilon = 1e-5;

	double S_N = 0, S_N_prev = 0, rel_chg = 0;
	int N = 0;

	double Nth_term = -1 / _x;    // sinTaylor3에서 sign_part와 pow_part와 fact_part를 모두 하나로 합쳐서 식을 작성해줌

	do {
		N++;
		S_N_prev = S_N;

		Nth_term *= -_x * _x / max((2 * N - 2) * (2 * N - 1), 1);
		S_N += Nth_term;

		rel_chg = fabs((S_N - S_N_prev) / S_N_prev);

	} while (N < N_max && rel_chg >= epsilon);

	return S_N;
}