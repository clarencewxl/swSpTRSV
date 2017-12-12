/*************************************************************************
	> File Name: my_solver.h
	> Author: Wang Xinliang
	> mail: clarencewxl@gmail.com
	> Created Time: Wed 22 Feb 2017 09:32:57 PM CST
 ************************************************************************/

#ifndef _MY_SOLVER_H_
#define _MY_SOLVER_H_

//#define PRODUCER_CONSUMER_ROWS      4
//#define PRODUCER_CONSUMER_ROWS_LOG2 2
#define PRODUCER_CONSUMER_ROWS_2NUM (PRODUCER_CONSUMER_ROWS-1)
#define CONSUMER_COLS               4
#define CONSUMER_COLS_LOG2          2
#define CONSUMER_COLS_2NUM          (CONSUMER_COLS-1)
#define PRODUCER_COLS               4
#define PRODUCER_COLS_LOG2          2
#define PRODUCER_COLS_2NUM          (PRODUCER_COLS-1)
#define CONSUMER_ROWS               PRODUCER_CONSUMER_ROWS
#define CONSUMER_ROWS_LOG2          PRODUCER_CONSUMER_ROWS_LOG2
#define CONSUMER_ROWS_2NUM          PRODUCER_CONSUMER_ROWS_2NUM
#define PRODUCER_ROWS               PRODUCER_CONSUMER_ROWS
#define PRODUCER_ROWS_LOG2          PRODUCER_CONSUMER_ROWS_LOG2
#define PRODUCER_ROWS_2NUM          PRODUCER_CONSUMER_ROWS_2NUM
#define CONSUMERS                   (CONSUMER_ROWS*CONSUMER_COLS)
#define CONSUMERS_LOG2              (CONSUMER_ROWS_LOG2+CONSUMER_COLS_LOG2)
#define CONSUMERS_2NUM              (CONSUMERS-1)
#define PRODUCERS                   (PRODUCER_ROWS*PRODUCER_COLS)
#define PRODUCERS_LOG2              (PRODUCER_ROWS_LOG2+PRODUCER_COLS_LOG2)
#define PRODUCERS_2NUM              (PRODUCERS-1)

#define VECTOR                    4
#define VECTOR_LOG2               2
#define VECTOR_2NUM               (VECTOR-1)

//#define CACHE_X_V4                1024
//#define CACHE_X_V4_LOG2           10
#define CACHE_X_V4_2NUM           (CACHE_X_V4-1)
#define CACHE_X                   (CACHE_X_V4*VECTOR)
#define CACHE_X_LOG2              (CACHE_X_V4_LOG2+VECTOR_LOG2)
#define CACHE_X_2NUM              (CACHE_X-1)

#define CONSUMER_CACHE            1024
#define CONSUMER_CACHE_LOG2       10
#define CONSUMER_CACHE_2NUM       (CONSUMER_CACHE-1)

#define CONSUMER_BUFFER_V4        256
#define CONSUMER_BUFFER_V4_LOG2   8
#define CONSUMER_BUFFER_V4_2NUM   (CONSUMER_BUFFER_V4-1)

#define PRODUCER_BUFFER_V4        256
#define PRODUCER_BUFFER_V4_LOG2   8
#define PRODUCER_BUFFER_V4_2NUM   (PRODUCER_BUFFER_V4-1)

#define PRODUCER_CACHE_LEVEL      1024
#define PRODUCER_CACHE_LEVEL_LOG2 10
#define PRODUCER_CACHE_LEVEL_2NUM (PRODUCER_CACHE_LEVEL-1)

#define PRODUCER_CACHE_IDX        64
#define PRODUCER_CACHE_IDX_LOG2   6
#define PRODUCER_CACHE_IDX_2NUM   (PRODUCER_CACHE_IDX-1)

#define PRODUCER_BUFFER_A         512
#define PRODUCER_BUFFER_A_LOG2    7
#define PRODUCER_BUFFER_A_2NUM    (PRODUCER_BUFFER_A-1)

#define MAX_SUB_LEVELS            (PRODUCER_ROWS+1)

#define BIG_BLOCK                 (PRODUCERS*CACHE_X_V4*VECTOR)
#define BIG_BLOCK_LOG2            (PRODUCERS_LOG2+CACHE_X_V4_LOG2+VECTOR_LOG2)
#define BIG_BLOCK_2NUM            (BIG_BLOCK-1)
#define BIG_BLOCK_V4              (PRODUCERS*CACHE_X_V4)
#define BIG_BLOCK_V4_LOG2         (PRODUCERS_LOG2+CACHE_X_V4_LOG2)
#define BIG_BLOCK_V4_2NUM         (BIG_BLOCK_V4-1)

typedef struct
{
	double *x,  *b;
	double *aa, *maa;
	int    *ai, *aj;
	int    *mai,*maj;
	int    *_idx, *idx[PRODUCERS], *level;
	int    *cnt[CONSUMERS];
	int    N, n, nnz, levels, mnnz;
} Matrix_Compress;

typedef struct
{
	double          *x,  *b;
	double          *aa, *maa;
	unsigned int    *ai, *aj;
	unsigned int    *mai,*maj;
	int             *_idx, *idx[PRODUCERS], *level, *plevel;
	int             *cnt[CONSUMERS];
	int             N, n, nnz, levels, plevels, mnnz;
    int             times;
} Matrix_Compress_General;

typedef struct
{
	int    n, level, nnz;
	int    *nz, *ai,  *aj;
	int    *ll, *ths, *ll_ths;
	int    *cnt[CONSUMERS];
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
	int    *cnt[CONSUMERS];
	int    N, n, nnz, levels, mnnz;
    int    times;
} Matrix_WXL;

#endif
