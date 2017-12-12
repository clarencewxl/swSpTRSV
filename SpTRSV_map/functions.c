/*************************************************************************
	> File Name: functions.c
	> Author: Wang Xinliang
	> mail: clarencewxl@gmail.com
	> Created Time: Mon 13 Feb 2017 06:35:07 PM CST
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <athread.h>
#include "my_slave.h"
#include "my_solver.h"

#define TIME(a,b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

extern SLAVE_FUN(stest)(void*);
extern SLAVE_FUN(Memcpy_MPE)(void*);
extern SLAVE_FUN(General_Solver_CSC_Multiple_MPE)(void*);

long long test()
{
	long long result = 0;
	athread_spawn(stest, &result);
	athread_join();
	return result;
}

int Memcpy(void *dst, void *src, int size)
{
	Memcpy_Info info;
	info.dst  = (char*)dst;
	info.src  = (char*)src;
	info.size = size;
	athread_spawn(Memcpy_MPE, &info);
	athread_join();
}

int Parallel_Solver_Compress_General(double *x, Matrix_Compress_General *m, double *b, int times)
{
	int i;
	int n = m->n;
	int N = m->N;
	int mnnz = m->mnnz;
	m->x = x;
	m->b = b;
    m->times = times;
    athread_spawn(General_Solver_CSC_Multiple_MPE, m);
	athread_join();
	return 0;
}

