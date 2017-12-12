/*************************************************************************
	> File Name: functions.h
	> Author: Wang Xinliang
	> mail: clarencewxl@gmail.com
	> Created Time: Mon 13 Feb 2017 06:40:47 PM CST
 ************************************************************************/

#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include "my_slave.h"
#include "my_solver.h"

int Parallel_Solver_Compress_General_R8B1024(double *x, Matrix_Compress_General *m, double *b, int times);
int Parallel_Solver_Compress_General_R4B512(double *x, Matrix_Compress_General *m, double *b, int times);
int Parallel_Solver_Compress_General_R2B64(double *x, Matrix_Compress_General *m, double *b, int times);
int Compress_General_R8B1024(Matrix_Compress_General *C, Matrix_Reorder *R);
int Compress_General_R4B512(Matrix_Compress_General *C, Matrix_Reorder *R);
int Compress_General_R2B64(Matrix_Compress_General *C, Matrix_Reorder *R);
#endif
