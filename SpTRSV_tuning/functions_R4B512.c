/*************************************************************************
	> File Name: functions.c
	> Author: Wang Xinliang
	> mail: clarencewxl@gmail.com
	> Created Time: Mon 13 Feb 2017 06:35:07 PM CST
 ************************************************************************/

#define PRODUCER_CONSUMER_ROWS      4
#define PRODUCER_CONSUMER_ROWS_LOG2 2
#define CACHE_X_V4                  512
#define CACHE_X_V4_LOG2             9 

#include <stdio.h>
#include <stdlib.h>
#include <athread.h>
#include "my_solver.h"
#include "my_slave.h"
#define TIME(a,b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

extern SLAVE_FUN(General_Solver_CSC_Multiple_MPE_R4B512)(void*);

int Compress_General_R4B512(Matrix_Compress_General *C, Matrix_Reorder *R)
{
	int i,j,k;
	int n   = C->n = R->n>>VECTOR_LOG2;
	int N   = C->N = R->n;
	int nnz = C->nnz = R->ai[N];
	int bn  = (n*VECTOR+BIG_BLOCK_2NUM)/BIG_BLOCK;
	C->mnnz = 0;	
	C->aa = malloc(nnz*sizeof(double));
	C->ai = malloc(nnz*sizeof(unsigned int));
	C->aj = malloc(nnz*sizeof(unsigned int));

    printf("autotuned PRODUCER_CONSUMER_ROWS:%d CACHE_X_V4:%d\n", PRODUCER_CONSUMER_ROWS, CACHE_X_V4);
	//printf("General Compress: The BN is %d\n", bn);
	
	int *tmp = malloc(n*sizeof(int));
	for(i = 0; i < n; i ++) tmp[i] = 0;
	for(i = 0; i < n*VECTOR; i ++)
		for(j = R->ai[i]; j < R->ai[i+1]; j ++)
			if(R->aj[j] < n*VECTOR)
			{
				if((i>>VECTOR_LOG2) != (R->aj[j]>>VECTOR_LOG2))
				{
					int level = tmp[i>>VECTOR_LOG2] + 1;
					if(level > tmp[R->aj[j]>>VECTOR_LOG2]) tmp[R->aj[j]>>VECTOR_LOG2] = level;
				}
			}
	int levels = 0;
	{
		int c = tmp[0];
		int t = 1;
		for(i = 1; i < n; i ++)
		{
			if(tmp[i] != c)
			{
				levels ++;
				t = 0;
			}
			else if(t == BIG_BLOCK_V4)
			{
				levels ++;
				t = 0;
			}
			else if((i%BIG_BLOCK_V4) == 0)
			{
				levels ++;
				t = 0;
			}
			c = tmp[i];
			t ++;
		}
		if(t > 0)
		{
			levels ++;
		}
	}
	//printf("General Compress: The BN Levels is %d\n", levels);
	C->levels = levels;
	C->level = malloc(levels*sizeof(int));
	//Set C->level
	{
		levels = 0;
		int c = tmp[0];
		int t = 1;
		for(i = 1; i < n; i ++)
		{
			if(tmp[i] != c)
			{
				C->level[levels++] = t;
				t = 0;
			}
			else if(t == BIG_BLOCK_V4)
			{
				C->level[levels++] = t;
				t = 0;
			}
			else if((i%BIG_BLOCK_V4) == 0)
			{
				C->level[levels++] = t;
				t = 0;
			}
			c = tmp[i];
			t ++;
		}
		if(t > 0)
		{
			C->level[levels++] = t;
		}
	}
	//for(i = 0; i < C->levels; i ++)
	//	printf("Level:%8d Values:%8d\n", i, C->level[i]);

    C->_cnt = malloc(CONSUMERS*levels*sizeof(int));
	//for(i = 0; i < CONSUMERS; i ++) C->cnt[i] = malloc(levels*sizeof(int));
	for(i = 0; i < CONSUMERS; i ++) C->cnt[i] = C->_cnt + i*levels;
	for(i = 0; i < CONSUMERS; i ++) for(j = 0; j < levels; j ++) C->cnt[i][j] = 0;
	
	int *len0 = malloc(     bn*bn *sizeof(int));
	int *len1 = malloc(PRODUCERS*PRODUCER_ROWS *bn*bn *sizeof(int));
	int *len2 = malloc(PRODUCERS*MAX_SUB_LEVELS*levels*sizeof(int));
	int *len3 = malloc(PRODUCERS  *bn*bn *sizeof(int));
	for(i = 0; i < bn*bn;       i ++) len0[i] = 0;
	for(i = 0; i < PRODUCERS*PRODUCER_ROWS *bn*bn;  i ++) len1[i] = 0;
	for(i = 0; i < PRODUCERS*MAX_SUB_LEVELS*levels; i ++) len2[i] = 0;
	for(i = 0; i < PRODUCERS*bn*bn;    i ++) len3[i] = 0;

	//Set len[012]
	int begin = 0, l = 0;
	for(l = 0; l < levels; l ++)
	{
		int nn = C->level[l];
		int offset = VECTOR*(begin+nn);
		for(i = 0; i < nn; i ++)
		{
			int pid = (i+begin)&PRODUCERS_2NUM;
			int id0 = pid >> PRODUCER_COLS_LOG2;
			for(j = (i+begin)*VECTOR; j < (i+1+begin)*VECTOR; j ++)
			{
				for(k = R->ai[j]; k < R->ai[j+1]; k ++)
				{
					int aj = R->aj[k];
					if(aj < n*VECTOR)
					{
						if((aj-offset) >= BIG_BLOCK)
						{
							len0[(aj>>BIG_BLOCK_LOG2)*bn+(j>>BIG_BLOCK_LOG2)] ++;
							int _pid = (j>>VECTOR_LOG2)&PRODUCERS_2NUM;
							int _id0 = _pid >> PRODUCER_COLS_LOG2;
							int _id1 = (aj>>(VECTOR_LOG2 + CONSUMER_COLS_LOG2))&CONSUMER_ROWS_2NUM;
							int _level = (_id1 + PRODUCER_ROWS - _id0) & PRODUCER_ROWS_2NUM;
							len1[((_pid+_level*PRODUCER_COLS)&PRODUCERS_2NUM)*bn*bn*PRODUCER_ROWS + \
                                (aj>>BIG_BLOCK_LOG2)*bn*PRODUCER_ROWS + \
                                (j>>BIG_BLOCK_LOG2)*PRODUCER_ROWS + \
                                _level] ++;
							len3[((aj>>VECTOR_LOG2)&CONSUMERS_2NUM)*bn*bn + (aj>>BIG_BLOCK_LOG2)*bn + (j>>BIG_BLOCK_LOG2)] ++;
						}
						else if(aj < (i+1+begin)*VECTOR)
						{
							len2[pid*MAX_SUB_LEVELS*levels+l*MAX_SUB_LEVELS+0] ++;
						}
						else
						{
							int id1 = (aj>>(VECTOR_LOG2 + CONSUMER_COLS_LOG2))&CONSUMER_ROWS_2NUM;
							int level = (id1 + PRODUCER_ROWS - id0) & PRODUCER_ROWS_2NUM;
							len2[((pid+level*PRODUCER_COLS)&PRODUCERS_2NUM)*MAX_SUB_LEVELS*levels+l*MAX_SUB_LEVELS+level+1] ++;
							C->cnt[(aj>>VECTOR_LOG2)&CONSUMERS_2NUM][l] ++;
						}
					}
					else
					{
						C->mnnz ++;
					}
				}
			}
		}
		begin += nn;
	}
	
    C->plevels = levels;
	for(i = 0; i < bn*bn; i ++) if(len0[i] > 0) C->plevels ++;
	C->plevel = malloc(C->plevels*sizeof(int));
	//Set C->plevel
	{
		int *plevel = C->plevel;
		int begin = 0;
		for(l = 0; l < levels; l ++)
		{
			int nn = C->level[l];
			*plevel++ = nn;
			begin += nn;
			if((begin&BIG_BLOCK_V4_2NUM) == 0 && l < (levels-1))
			{
				int *ll = len0+(begin>>BIG_BLOCK_V4_LOG2)*bn;
				for(i = 0; i < bn; i ++)
					if(ll[i] > 0) *plevel++ = -i;
			}
		}
	}
	//Check C->plevel
	//{
	//	for(i = 0; i < C->plevels; i ++)
	//		printf("Compress General check plevel %4d %8d\n", i, C->plevel[i]);
	//}

	int plevels = C->plevels;
	C->_idx = malloc((PRODUCERS*MAX_SUB_LEVELS*plevels+1)*sizeof(int)); //C->_idx[32][plevels][9]
	for(i = 0; i <  PRODUCERS;                        i ++) C->idx[i]  = C->_idx + i*MAX_SUB_LEVELS*plevels;
	for(i = 0; i <= PRODUCERS*MAX_SUB_LEVELS*plevels; i ++) C->_idx[i] = 0;

	unsigned int **ai1 = malloc(PRODUCERS*PRODUCER_ROWS*  bn*bn*sizeof(unsigned int*));     //[32][bn][bn][8]
	unsigned int **ai2 = malloc(PRODUCERS*MAX_SUB_LEVELS*levels*sizeof(unsigned int*));    //[32][levels][9]
	unsigned int **aj1 = malloc(PRODUCERS*PRODUCER_ROWS*  bn*bn*sizeof(unsigned int*));     //[32][bn][bn][8]
	unsigned int **aj2 = malloc(PRODUCERS*MAX_SUB_LEVELS*levels*sizeof(unsigned int*));    //[32][levels][9]
	double       **aa1 = malloc(PRODUCERS*PRODUCER_ROWS*bn*bn*sizeof(double*));  //[32][bn][bn][8]
	double       **aa2 = malloc(PRODUCERS*MAX_SUB_LEVELS*levels*sizeof(double*)); //[32][levels][9]
    
    //printf("Compree General check malloc: (%lld) %p %p %p %p %p %p\n", PRODUCERS*PRODUCER_ROWS*  bn*bn*sizeof(unsigned int*), ai1, ai2, aj1, aj2, aa1, aa2);
	//Set a[ija][12] and idx
	{
		int   *_idx = C->_idx + 1;
		unsigned int    *_ai = C->ai;
		unsigned int    *_aj = C->aj;
		double *_aa = C->aa;
		for(i = 0; i < PRODUCERS; i ++)
		{
			int begin = 0;
			for(l = 0; l < levels; l ++)
			{
				int nn = C->level[l];
				for(j = 0; j < MAX_SUB_LEVELS; j ++)
				{
					ai2[i*levels*MAX_SUB_LEVELS+l*MAX_SUB_LEVELS+j] = _ai; _ai += len2[i*levels*MAX_SUB_LEVELS+l*MAX_SUB_LEVELS+j];
					aj2[i*levels*MAX_SUB_LEVELS+l*MAX_SUB_LEVELS+j] = _aj; _aj += len2[i*levels*MAX_SUB_LEVELS+l*MAX_SUB_LEVELS+j];
					aa2[i*levels*MAX_SUB_LEVELS+l*MAX_SUB_LEVELS+j] = _aa; _aa += len2[i*levels*MAX_SUB_LEVELS+l*MAX_SUB_LEVELS+j];
					*_idx++ = len2[i*levels*MAX_SUB_LEVELS+l*MAX_SUB_LEVELS+j];
				}
				begin += nn;
				if((begin&BIG_BLOCK_V4_2NUM) == 0 && l < (levels-1))
				{
					int ii = begin>>BIG_BLOCK_V4_LOG2;
					for(k = 0; k < bn; k ++)
					{
						for(j = 0; j < PRODUCER_ROWS; j ++)
						{
							ai1[i*bn*bn*PRODUCER_ROWS + ii*bn*PRODUCER_ROWS + k*PRODUCER_ROWS + j] = _ai; _ai += len1[i*bn*bn*PRODUCER_ROWS + ii*bn*PRODUCER_ROWS + k*PRODUCER_ROWS + j];
							aj1[i*bn*bn*PRODUCER_ROWS + ii*bn*PRODUCER_ROWS + k*PRODUCER_ROWS + j] = _aj; _aj += len1[i*bn*bn*PRODUCER_ROWS + ii*bn*PRODUCER_ROWS + k*PRODUCER_ROWS + j];
							aa1[i*bn*bn*PRODUCER_ROWS + ii*bn*PRODUCER_ROWS + k*PRODUCER_ROWS + j] = _aa; _aa += len1[i*bn*bn*PRODUCER_ROWS + ii*bn*PRODUCER_ROWS + k*PRODUCER_ROWS + j];
						}
						if(len0[ii*bn+k] > 0)
						{
							*_idx++ = 0;
							for(j = 0; j < PRODUCER_ROWS; j ++)
								*_idx++ = len1[i*bn*bn*PRODUCER_ROWS + ii*bn*PRODUCER_ROWS + k*PRODUCER_ROWS + j];
						}
						C->cnt[i][l] += len3[i*bn*bn+ii*bn+k];
					}

				}
			}
		}
		for(i = 0; i < PRODUCERS*MAX_SUB_LEVELS*plevels; i ++) C->_idx[i+1] += C->_idx[i];
	}

	C->mnnz += (R->ai[N]-R->ai[n*VECTOR]);
	//for(i = 0; i < PRODUCERS; i ++) 
    //    printf("Compress Check idx: Pid %2d From %8d  To %8d (%8d)\n", i, C->idx[i][0], C->idx[i][MAX_SUB_LEVELS*plevels], C->idx[i][MAX_SUB_LEVELS*plevels]-C->idx[i][0]); 

	//Init m information
	C->maa = NULL;
	C->mai = NULL;
	C->maj = NULL;
	if(C->mnnz > 0)
	{
		C->maa = malloc(C->mnnz*sizeof(double));
		C->mai = malloc(C->mnnz*sizeof(unsigned int));
		C->maj = malloc(C->mnnz*sizeof(unsigned int));
	}
	//printf("Compress Check nnz and last idx: %8d, %8d, %8d\n", C->nnz, C->_idx[PRODUCERS*MAX_SUB_LEVELS*plevels], C->mnnz);
	double       *maa = C->maa;
	unsigned int *mai = C->mai;
	unsigned int *maj = C->maj;
	begin = 0, l = 0;
	for(l = 0; l < levels; l ++)
	{
		int nn = C->level[l];
		int base = (begin+nn)*VECTOR;
		int offset = VECTOR*(begin+nn);
		for(i = 0; i < nn; i ++)
		{
			int pid = (i+begin)&PRODUCERS_2NUM;
			int id0 = pid>>PRODUCER_COLS_LOG2;
			for(j = (i+begin)*VECTOR; j < (i+1+begin)*VECTOR; j ++)
			{
				for(k = R->ai[j]; k < R->ai[j+1]; k ++)
				{
					int ai = j;
					int aj = R->aj[k];
					double       aa = R->a [k];
					if(aj < n*VECTOR)
					{
						if((aj-offset) >= BIG_BLOCK)
						{
							int _pid = (j>>VECTOR_LOG2)&PRODUCERS_2NUM;
							int _id0 = _pid >> PRODUCER_COLS_LOG2;
							int _id1 = (aj>>(VECTOR_LOG2+CONSUMER_COLS_LOG2))&CONSUMER_ROWS_2NUM;
							unsigned int _level = (_id1+PRODUCER_ROWS-_id0)&PRODUCER_ROWS_2NUM;
							unsigned int aai = (((ai&BIG_BLOCK_2NUM)>>(PRODUCERS_LOG2+VECTOR_LOG2))<<VECTOR_LOG2)+(ai&VECTOR_2NUM);
							unsigned int aaj =    aj&BIG_BLOCK_2NUM;
							*aj1[((_pid+_level*PRODUCER_COLS)&PRODUCERS_2NUM)*bn*bn*PRODUCER_ROWS + (aj>>BIG_BLOCK_LOG2)*bn*PRODUCER_ROWS + (j>>BIG_BLOCK_LOG2)*PRODUCER_ROWS + _level]++ = aaj + (aai<<BIG_BLOCK_LOG2) + (_level<<(BIG_BLOCK_LOG2+CACHE_X_LOG2));
							//*aj1[((_pid+_level*PRODUCER_COLS)&PRODUCERS_2NUM)*bn*bn*PRODUCER_ROWS + (aj>>BIG_BLOCK_LOG2)*bn*PRODUCER_ROWS + (j>>BIG_BLOCK_LOG2)*PRODUCER_ROWS + _level]++ = aaj + (aai<<BIG_BLOCK_LOG2);
							*ai1[((_pid+_level*PRODUCER_COLS)&PRODUCERS_2NUM)*bn*bn*PRODUCER_ROWS + (aj>>BIG_BLOCK_LOG2)*bn*PRODUCER_ROWS + (j>>BIG_BLOCK_LOG2)*PRODUCER_ROWS + _level]++ = ai;
							*aa1[((_pid+_level*PRODUCER_COLS)&PRODUCERS_2NUM)*bn*bn*PRODUCER_ROWS + (aj>>BIG_BLOCK_LOG2)*bn*PRODUCER_ROWS + (j>>BIG_BLOCK_LOG2)*PRODUCER_ROWS + _level]++ = aa;
						}
						else if(aj < (i+1+begin)*VECTOR)
						{
							unsigned int aai = ((((ai>>VECTOR_LOG2)-begin)>>PRODUCERS_LOG2)<<VECTOR_LOG2)+(ai&VECTOR_2NUM);
							unsigned int aaj = ((((aj>>VECTOR_LOG2)-begin)>>PRODUCERS_LOG2)<<VECTOR_LOG2)+(aj&VECTOR_2NUM);
							*aj2[pid*MAX_SUB_LEVELS*levels+l*MAX_SUB_LEVELS+0]++ = (aai<<BIG_BLOCK_LOG2)+aaj;
							*ai2[pid*MAX_SUB_LEVELS*levels+l*MAX_SUB_LEVELS+0]++ = ai;
							*aa2[pid*MAX_SUB_LEVELS*levels+l*MAX_SUB_LEVELS+0]++ = aa;
							//if(ai == aj) *aa2[pid*MAX_SUB_LEVELS*levels+l*MAX_SUB_LEVELS+0]++ = 1-aa;
							//else         *aa2[pid*MAX_SUB_LEVELS*levels+l*MAX_SUB_LEVELS+0]++ = aa;
						}
						else
						{
							int id1 = (aj>>(VECTOR_LOG2+CONSUMER_COLS_LOG2))&CONSUMER_ROWS_2NUM;
							unsigned int level = (id1+PRODUCER_ROWS-id0)&PRODUCER_ROWS_2NUM;
							unsigned int aai = ((((ai>>VECTOR_LOG2)-begin)>>PRODUCERS_LOG2)<<VECTOR_LOG2)+(ai&VECTOR_2NUM);
							unsigned int aaj = aj-base;
							*aj2[((pid+level*PRODUCER_COLS)&PRODUCERS_2NUM)*MAX_SUB_LEVELS*levels+l*MAX_SUB_LEVELS+level+1]++ = aaj + (aai<<BIG_BLOCK_LOG2) + (level<<(BIG_BLOCK_LOG2+CACHE_X_LOG2));
							//*aj2[((pid+level*PRODUCER_COLS)&PRODUCERS_2NUM)*MAX_SUB_LEVELS*levels+l*MAX_SUB_LEVELS+level+1]++ = aaj + (aai<<BIG_BLOCK_LOG2);
							*ai2[((pid+level*PRODUCER_COLS)&PRODUCERS_2NUM)*MAX_SUB_LEVELS*levels+l*MAX_SUB_LEVELS+level+1]++ = ai;
							*aa2[((pid+level*PRODUCER_COLS)&PRODUCERS_2NUM)*MAX_SUB_LEVELS*levels+l*MAX_SUB_LEVELS+level+1]++ = aa;
						}
					}
					else
					{
						int _id = (ai>>VECTOR_LOG2)&PRODUCERS_2NUM;
						int _offset = _id*(n>>PRODUCERS_LOG2) + (_id < (n&PRODUCERS_2NUM) ? _id : (n&PRODUCERS_2NUM)) + (ai>>(PRODUCERS_LOG2+VECTOR_LOG2));
						*mai ++ = _offset*VECTOR+(ai&VECTOR_2NUM);
						*maj ++ = aj;
						*maa ++ = aa;
					}
				}
			}	
		}
		begin += nn;
	}
	for(i = VECTOR*n; i < N; i ++)
	{
		for(j = R->ai[i]; j < R->ai[i+1]; j ++)
		{
			if(R->aj[j] == i)
			{
				*mai ++ = i;
				*maj ++ = i;
				*maa ++ = R->a[R->ai[i]];
				//*maa ++ = 1-R->a[R->ai[i]];
			}
			else
			{
				*mai ++ = i;
				*maj ++ = R->aj[j];
				*maa ++ = R->a [j];
			}
		}
	}
	//Check a[ija][12]
	{
		unsigned int    *_ai = C->ai;
		unsigned int    *_aj = C->aj;
		double *_aa = C->aa;
		for(i = 0; i < PRODUCERS; i ++)
		{
			int begin = 0;
			for(l = 0; l < levels; l ++)
			{
				int nn = C->level[l];
				for(j = 0; j < MAX_SUB_LEVELS; j ++)
				{
					_ai += len2[i*levels*MAX_SUB_LEVELS+l*MAX_SUB_LEVELS+j];
					_aj += len2[i*levels*MAX_SUB_LEVELS+l*MAX_SUB_LEVELS+j];
					_aa += len2[i*levels*MAX_SUB_LEVELS+l*MAX_SUB_LEVELS+j];
					if(ai2[i*levels*MAX_SUB_LEVELS+l*MAX_SUB_LEVELS+j] != _ai) printf("Compress General Error ai %2d %4d %1d\n", i, l, j);
					if(aj2[i*levels*MAX_SUB_LEVELS+l*MAX_SUB_LEVELS+j] != _aj) printf("Compress General Error aj %2d %4d %1d\n", i, l, j);
					if(aa2[i*levels*MAX_SUB_LEVELS+l*MAX_SUB_LEVELS+j] != _aa) printf("Compress General Error aa %2d %4d %1d\n", i, l, j);
				}
				begin += nn;
				if((begin&BIG_BLOCK_V4_2NUM) == 0 && l < (levels-1))
				{
					int ii = begin>>BIG_BLOCK_V4_LOG2;
					for(k = 0; k < bn; k ++)
					{
						for(j = 0; j < PRODUCER_ROWS; j ++)
						{
							_ai += len1[i*bn*bn*PRODUCER_ROWS + ii*bn*PRODUCER_ROWS + k*PRODUCER_ROWS + j];
							_aj += len1[i*bn*bn*PRODUCER_ROWS + ii*bn*PRODUCER_ROWS + k*PRODUCER_ROWS + j];
							_aa += len1[i*bn*bn*PRODUCER_ROWS + ii*bn*PRODUCER_ROWS + k*PRODUCER_ROWS + j];
							if(ai1[i*bn*bn*PRODUCER_ROWS + ii*bn*PRODUCER_ROWS + k*PRODUCER_ROWS + j] != _ai) printf("Compress General Error ai %2d %4d %8d %1d\n", i, l, k, j);
							if(aj1[i*bn*bn*PRODUCER_ROWS + ii*bn*PRODUCER_ROWS + k*PRODUCER_ROWS + j] != _aj) printf("Compress General Error aj %2d %4d %8d %1d\n", i, l, k, j);
							if(aa1[i*bn*bn*PRODUCER_ROWS + ii*bn*PRODUCER_ROWS + k*PRODUCER_ROWS + j] != _aa) printf("Compress General Error aa %2d %4d %8d %1d\n", i, l, k, j);
						}
					}
				}
			}
		}
	}
	free(len0);free(len1);free(len2);free(len3);
	free(ai1); free(aj1); free(aa1);
	free(ai2); free(aj2); free(aa2);
	free(tmp);
	return 0;
}

int Parallel_Solver_Compress_General_R4B512(double *x, Matrix_Compress_General *m, double *b, int times)
{
	int i;
	int n = m->n;
	int N = m->N;
	int mnnz = m->mnnz;
	m->x = x;
	m->b = b;
    m->times = times;
    athread_spawn(General_Solver_CSC_Multiple_MPE_R4B512, m);
	athread_join();
	return 0;
}

