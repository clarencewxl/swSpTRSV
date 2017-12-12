/*************************************************************************
	> File Name: my_solver.h
	> Author: Wang Xinliang
	> mail: clarencewxl@gmail.com
	> Created Time: Wed 22 Feb 2017 09:32:57 PM CST
 ************************************************************************/

#ifndef _MY_SOLVER_H_
#define _MY_SOLVER_H_


typedef struct
{
	double *x,  *b;
	double *aa, *maa;
	int    *ai, *aj;
	int    *mai,*maj;
	int    *_idx, *idx[32], *level;
	int    *cnt[32];
	int    N, n, nnz, levels, mnnz;
} Matrix_Compress;

typedef struct
{
	double          *x,  *b;
	double          *aa, *maa;
	unsigned int    *ai, *aj;
	unsigned int    *mai,*maj;
	int             *_idx, *idx[32], *level, *plevel;
	int             *_cnt, *cnt[32];
	int             N, n, nnz, levels, plevels, mnnz;
    int             times;
} Matrix_Compress_General;

typedef struct
{
	int    n, level, nnz;
	int    *nz, *ai,  *aj;
	int    *ll, *ths, *ll_ths;
	int    *cnt[32];
	double *aa;
	double *x;
	double *b;
} Matrix_Special;

typedef struct
{
	double *x,  *b;
	double *aa, *maa;
	int    *ai, *aj;
	int    *mai, *maj;
	int    *idx, *level;
	int    *cnt[32];
	int    N, n, nnz, levels, mnnz;
} Matrix_WXL;

typedef struct
{
	double *a;
	int *ai;
	int *aj;
	int *nz;
	int n;
} Matrix;

typedef struct
{
	double *a;
	int *ai;
	int *aj;
	int *levels;
    int *m2r;
	int *r2m;
	int n, level;
} Matrix_Reorder;

typedef struct
{
    double *x, *b;
    double *aa;
    int *ai, *aj;
    int *levels;
    int N, L;
    int times;
} Matrix_LS;

#endif
