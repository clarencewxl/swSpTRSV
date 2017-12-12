/*************************************************************************
	> File Name: functions.h
	> Author: Wang Xinliang
	> mail: clarencewxl@gmail.com
	> Created Time: Mon 13 Feb 2017 06:40:47 PM CST
 ************************************************************************/

#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include "my_solver.h"

long long test();
int Init_Solver(int *ai, int *aj, int N);
int Finalize_Solver();
int Solver_CSC_Single  (double *x, double *a, int *ai, int *aj, int *nz, double *b, int N);
int Solver_CSC_Multiple(double *x, Matrix_Special *s, double *b);
int Parallel_Solver_WXL(double *x, Matrix_WXL *m, double *b);
int Parallel_Solver_Compress_General(double *x, Matrix_Compress_General *m, double *b, int times);

#endif
