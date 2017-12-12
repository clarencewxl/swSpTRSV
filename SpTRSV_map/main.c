/*************************************************************************
	> File Name: main.c
	> Author: Wang Xinliang
	> mail: clarencewxl@gmail.com
	> Created Time: Mon 13 Feb 2017 06:35:07 PM CST
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <athread.h>
#include <math.h>
#include "functions.h"
#define TIME(a,b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

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
	int *m2r;
	int *r2m;
	int n, level;
} Matrix_Reorder;

int _cmp(const void *a , const void *b )
{
	return *(int *)a - *(int *)b;
}
int _cmp_levels(const void *_a , const void *_b )
{
	int *a = (int*)_a;
	int *b = (int*)_b;
	int result = a[0]-b[0];
	if(result != 0) return result;
	else			return a[1]-b[1];
}

int Compress(Matrix_Compress *C, Matrix_Reorder *R)
{
	int i,j,k;
	int n   = C->n = R->n>>2;
	int N   = C->N = R->n;
	int nnz = C->nnz = R->ai[N];
	C->mnnz = 0;	
	C->aa = malloc(nnz*sizeof(double));
	C->ai = malloc(nnz*sizeof(int));
	C->aj = malloc(nnz*sizeof(int));
	
	int *tmp = malloc(n*sizeof(int));
	for(i = 0; i < n; i ++) tmp[i] = 0;
	for(i = 0; i < n*4; i ++)
		for(j = R->ai[i]; j < R->ai[i+1]; j ++)
			if(R->aj[j] < n*4)
			{
				if((i>>2) != (R->aj[j]>>2))
				{
					int level = tmp[i>>2] + 1;
					if(level > tmp[R->aj[j]>>2]) tmp[R->aj[j]>>2] = level;
				}
			}
	//for(i = 0; i < n; i ++) printf("Check level %8d(%8d-%8d) %8d\n", i, 4*i, 4*i+3, tmp[i]);
	int levels = C->levels = tmp[n-1] + 1;
	printf("Compress Check: Total levels %d\n", levels);
	C->level = malloc(levels*sizeof(int));
	C->_idx  = malloc((32*9*levels+1)*sizeof(int));
	for(i = 0; i < 32; i ++) C->idx[i] = C->_idx + i*9*levels;
	for(i = 0; i < 32; i ++) C->cnt[i] = malloc(levels*sizeof(int));
	for(i = 0; i < 32; i ++)
		for(j = 0; j < levels; j ++) C->cnt[i][j] = 0;
	for(i = 0; i <  levels; i ++) C->level[i] = 0;
	for(i = 0; i <= 32*9*levels; i ++) C->_idx[i] = 0;
	//Set C->level
	{
		int *l = C->level-1;
		int f = -1;
		for(i = 0; i < n; i ++)
		{
			if(tmp[i] != f)
			{
				f = tmp[i];
				l ++;
			}
			l[0] ++;
		}

	}
	for(i = 0; i < levels; i ++)
	{
		printf("Compress Check level %8d has %8d values\n", i, C->level[i]);
	}
	for(i = 0; i < levels; i ++)
	{
		if(C->level[i] > 32768)
		{
			printf("Compress ERROR 1: Level %d has to many blocks (%d)\n", i, C->level[i]);
			exit(-1);
		}
	}

	//Set C->idx
	int begin = 0, l = 0;
	int signal = 0;
	for(l = 0; l < levels; l ++)
	{
		int nn = C->level[l];
		//int dnn = nn & 0x1F;
		//int pnn = nn >> 5;
		int offset = 4*(begin+nn);
		int *idx[32];
		for(i = 0; i < 32; i ++) idx[i] = C->idx[i] + 1 + l*9;
		for(i = 0; i < nn; i ++)
		{
			int pid = (i+begin)&0x1F;
			int id0 = pid>>2;
			for(j = (i+begin)*4; j < (i+1+begin)*4; j ++)
			{
				for(k = R->ai[j]; k < R->ai[j+1]; k ++)
				{
					int aj = R->aj[k];
					if(aj < n*4)
					{
						if((aj-offset) >= 131072)
						{
							signal ++;
						}
						if(aj < (i+1+begin)*4)
						{
							idx[pid][0] ++;
						}
						else
						{
							int id1 = (aj>>4)&0x07;
							int level = (id1+8-id0)&0x07;
							idx[(pid+level*4)&0x1F][level+1]++;
							C->cnt[(aj>>2)&0x1F][l] ++;
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
	if(signal > 0)
	{
		printf("Compress ERROR 2: signal:%8d nnz:%8d ratio:%.6lf%\n", signal, nnz, ((double)signal)/((double)nnz)*100.0);
		exit(-1);
	}
	C->mnnz += (R->ai[N]-R->ai[n*4]);
	for(i = 0; i < 32*9*levels; i ++) C->_idx[i+1] += C->_idx[i];
	for(i = 0; i < 32; i ++) printf("Compress Check idx: Pid %2d From %8d  To %8d (%8d)\n", i, C->idx[i][0], C->idx[i][9*levels], C->idx[i][9*levels]-C->idx[i][0]); 

	//Init m information
	C->maa = NULL;
	C->mai = NULL;
	C->maj = NULL;
	if(C->mnnz > 0)
	{
		C->maa = malloc(C->mnnz*sizeof(double));
		C->mai = malloc(C->mnnz*sizeof(int));
		C->maj = malloc(C->mnnz*sizeof(int));
	}
	printf("Compress Check nnz and last idx: %8d, %8d, %8d\n", C->nnz, C->_idx[32*9*levels], C->mnnz);

	//Set C->aa, C->ai, C->aj
	double *maa = C->maa;
	int    *mai = C->mai;
	int    *maj = C->maj;
	begin = 0, l = 0;
	for(l = 0; l < levels; l ++)
	{
		int nn = C->level[l];
		int base = (begin+nn)*4;
		int offset = 4*(begin+nn);
		int    *lai[32][9];
		int    *laj[32][9];
		double *laa[32][9];
		for(i = 0; i < 32; i ++)
			for(j = 0; j < 9; j ++)
			{
				lai[i][j] = C->ai + C->idx[i][l*9+j];
				laj[i][j] = C->aj + C->idx[i][l*9+j];
				laa[i][j] = C->aa + C->idx[i][l*9+j];
			}
		for(i = 0; i < nn; i ++)
		{
			int pid = (i+begin)&0x1F;
			int id0 = pid>>2;
			for(j = (i+begin)*4; j < (i+1+begin)*4; j ++)
			{
				for(k = R->ai[j]; k < R->ai[j+1]; k ++)
				{
					int    ai = j;
					int    aj = R->aj[k];
					double aa = R->a [k];
					if(aj < (i+1+begin)*4)
					{
						int aai = ((((ai>>2)-begin)>>5)<<2)+(ai&0x03);
						int aaj = ((((aj>>2)-begin)>>5)<<2)+(aj&0x03);
						*laj[pid][0]++ = (aai<<17)+aaj;
						//*lai[pid][0]++ = ai;
						//*laj[pid][0]++ = aj;
						if(ai == aj) *laa[pid][0]++ = 1-aa;
						else         *laa[pid][0]++ = aa;
					}
					else if(aj < 4*n)
					{
						int id1 = (aj>>4)&0x07;
						int level = (id1+8-id0)&0x07;
						int aai = ((((ai>>2)-begin)>>5)<<2)+(ai&0x03);
						int aaj = aj-base;
						*laj[(pid+level*4)&0x1F][level+1]++ = (aai<<17)+aaj;
						//*lai[(pid+level*4)&0x1F][level+1]++ = ((((ai>>2)-begin)>>5)<<2)+(ai&0x03);
						//*laj[(pid+level*4)&0x1F][level+1]++ = aj;
						*laa[(pid+level*4)&0x1F][level+1]++ = aa;
					}
					else
					{
						int _id = (ai>>2)&0x1F;
						int _offset = _id*(n>>5) + (_id < (n&0x1F) ? _id : (n&0x1F)) + (ai>>7);
						*mai ++ = _offset*4+(ai&0x03);
						*maj ++ = aj;
						*maa ++ = aa;
					}
					
				}
			}	
		}
		for(i = 0; i < 32; i ++)
			for(j = 0; j < 9; j ++)
			{
				//if(lai[i][j] != C->ai + C->idx[i][l*9+j+1]) printf("Error coding in ai %d %d %d\n", l, i, j);
				if(laj[i][j] != C->aj + C->idx[i][l*9+j+1]) printf("Error coding in aj %d %d %d\n", l, i, j);
				if(laa[i][j] != C->aa + C->idx[i][l*9+j+1]) printf("Error coding in aa %d %d %d\n", l, i, j);
			}
		begin += nn;
	}
	for(i = 4*n; i < N; i ++)
	{
		for(j = R->ai[i]; j < R->ai[i+1]; j ++)
		{
			if(R->aj[j] == i)
			{
		*mai ++ = i;
		*maj ++ = i;
		*maa ++ = 1-R->a[R->ai[i]];
			}
			else
			{
			*mai++ = i;
			*maj++ = R->aj[j];
			*maa++ = R->a [j];
			}
		}

	}

	free(tmp);
	return 0;
}

int Compress_General(Matrix_Compress_General *C, Matrix_Reorder *R)
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

	for(i = 0; i < CONSUMERS; i ++) C->cnt[i] = malloc(levels*sizeof(int));
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

double Reorder(Matrix_Reorder *R, Matrix *M)
{
	int i,j,k;
	int n   = R->n = M->n;
	int nnz = M->ai[n] - M->ai[0];
	int max_level = 0;
	int *levels = (int*)malloc(2*n*sizeof(int));
	R->m2r   = (int*)malloc(n*sizeof(int));
	R->r2m   = (int*)malloc(n*sizeof(int));
	R->ai    = (int*)   malloc((n+1)* sizeof(int));
	R->aj    = (int*)   malloc( nnz * sizeof(int));
	R->a     = (double*)malloc( nnz * sizeof(double));
	for(i = 0; i < n; i ++)
	{
		levels[2*i+1] = i;
		levels[2*i+0] = 0;
	}
	for(i = 0; i < n; i ++)
	{
		for(j = M->ai[i]; j < M->ai[i+1]; j ++)
		{
			if(i != M->aj[j])
			{
				int level = levels[2*i] + 1;
				if(level > levels[2*M->aj[j]])
					levels[2*M->aj[j]] = level;
				if(level > max_level) max_level = level;
			}
		}
	}
	max_level ++;
	qsort(levels, n, sizeof(int)*2, _cmp_levels);
	for(i = 0; i < n; i ++)
	{
		R->m2r[levels[2*i+1]] = i;
		R->r2m[i] = levels[2*i+1];
	}
	
	int *cnt    = malloc((max_level+1)*sizeof(int));
	for(i = 0; i <= max_level; i ++) cnt[i] = 0;
	for(i = 0; i < n; i ++) cnt[levels[2*i]+1] ++;
	for(i = 0; i < max_level; i ++) cnt[i+1] += cnt[i];
	
	for(i = 0; i < n; i ++) levels[2*i] = 0;

	for(i = 0; i < max_level; i ++)
	{
		for(j = cnt[i]; j < cnt[i+1]; j ++)
		{
			int jj = levels[2*j+1];
			for(k = M->ai[jj]; k < M->ai[jj+1]; k ++)
			{
				if(levels[2*R->m2r[M->aj[k]]] < j)
					levels[2*R->m2r[M->aj[k]]] = j;			
			}
		}
		if(i < (max_level-1))
		{
			qsort(levels+2*cnt[i+1], cnt[i+2]-cnt[i+1], sizeof(int)*2, _cmp_levels);
		}
	}
	for(i = 0; i < n; i ++)
	{
		R->m2r[levels[2*i+1]] = i;
		R->r2m[i] = levels[2*i+1];
	}

	free(cnt);


	R->ai[0] = 0;
	for(i = 0; i < n; i ++)
	{
		int m  = R->r2m[i];
		int nz = M->ai[m+1]-M->ai[m];
		R->ai[i+1] = R->ai[i] + nz;
		for(j = 0; j < nz; j ++)
		{
			R->aj[R->ai[i]+j] = R->m2r[M->aj[M->ai[m]+j]];
			R->a [R->ai[i]+j] = M->a[M->ai[m]+j];
		}
	}
	//for(i = 0; i < n; i ++)
	//{
	//	for(j = R->ai[i]+1; j < R->ai[i+1]-1; j ++)
	//	{
	//		for(k = 0; k < R->ai[i+1]-j-1; k ++)
	//		{
	//			if(R->aj[j+k] > R->aj[j+1+k])
	//			{
	//				int t = R->aj[j+k];
	//				R->aj[j+k] = R->aj[j+1+k];
	//				R->aj[j+1+k] = t;
	//				double v = R->a[j+k];
	//				R->a[j+k] = R->a[j+1+k];
	//				R->a[j+1+k] = v;
	//			}
	//		}
	//	}
	//}
	for(i = 0; i < n; i ++)
	{
		levels[2*i+1] = i;
		levels[2*i+0] = 0;
	}
	for(i = 0; i < n; i ++)
	{
		for(j = R->ai[i]+1; j < R->ai[i+1]; j ++)
		{
			int level = levels[2*i] + 1;
			if(level > levels[2*R->aj[j]])
				levels[2*R->aj[j]] = level;
		}
	}
	for(i = 0; i < n-1; i ++)
	{
		if(levels[2*i+2] < levels[2*i])
		{
			printf("Warning %d\n", i);
		}
	}

	free(levels);
#if 0
	//Check Reorder
	{
		double *data = malloc(4*n*sizeof(double));
		double *x0 = data + 0*n;
		double *x1 = data + 1*n;
		double *b0 = data + 2*n;
		double *b1 = data + 3*n;
		int    *m2r = R->m2r;
		for(i = 0; i < n; i ++)
		{
			x1[m2r[i]] = x0[i] = 0.01*i;
			b0[i] = b1[i] = 0;
		}
		printf("Begin to check reorder\n");
		for(i = 0; i < n; i ++)
		{
			for(j = M->ai[i]; j < M->ai[i+1]; j ++)
			{
				b0[M->aj[j]] += M->a[j] * x0[M->aj[j]];
			}
		}
		for(i = 0; i < n; i ++)
		{
			for(j = R->ai[i]; j < R->ai[i+1]; j ++)
			{
				b1[R->aj[j]] += R->a[j] * x1[R->aj[j]];
			}	
		}
		for(i = 0; i < n; i ++)
		{
			double tmp = b0[i]-b1[m2r[i]];
			if(fabs(tmp) > 1e-10)
			{
				printf("Wrong: The Reorder method is wrong in M:%d and R:%d (%.14lf %.14lf)\n", i, m2r[i], b0[i], b1[m2r[i]]);
				break;
			}

		}
		printf("Check reorder is over\n");
		free(data);
	}
#endif
	return (1.0*nnz)/(1.0*max_level);;
}


void check(double *A, double *B, int N)
{
	int i;
	for(i = 0; i < N; i ++)
	{
		if(fabs(A[i]-B[i])/MAX(fabs(A[i]),1.0) > 1e-12)
		{
			printf("Warning: The first wrong %d(%.16lf %.16lf)\n", i, A[i], B[i]);
			exit(-1);
		}
	}
}

void special_check(double *A, double *B, int N)
{
	int i;
	int n = N>>VECTOR_LOG2;
	int pn = n>>PRODUCERS_LOG2;
	int dn = n&PRODUCERS_2NUM;
	for(i = 0; i < n*VECTOR; i ++)
	{
		int id = (i>>VECTOR_LOG2)&PRODUCERS_2NUM;
		int offset = id*pn + (id < dn ? id : dn) + (i>>(VECTOR_LOG2+PRODUCERS_LOG2));
		int j = offset*VECTOR+(i&VECTOR_2NUM);
		if(fabs(A[i]-B[j])/MAX(fabs(A[i]),1.0) > 1e-12)
		{
			printf("Warning: The first wrong %d(%.16lf %.16lf)\n", i, A[i], B[j]);
			exit(-1);
		}
	}
	for(i; i < N; i ++)
	{
		int j = i;
		if(fabs(A[i]-B[j])/MAX(fabs(A[i]),1.0) > 1e-12)
		{
			printf("Warning: The first wrong %d(%.16lf %.16lf)\n", i, A[i], B[j]);
			exit(-1);
		}
	}
}

int Init_Random(Matrix *m, int n, int max)
{
	int i,j,k;
	int *signal = malloc(n*sizeof(int));
	m->n = n;
	m->ai = malloc((n+1)*sizeof(int));
	//Init ai
	m->ai[0] = 0;
	for(i = 1; i <= n; i ++)
		m->ai[i] = m->ai[i-1] + (rand()%MIN(max+1,n-i+1)) + 1;
	//Init aj
	for(i = 0; i < n; i ++) signal[i] = 1;
	m->aj = malloc(m->ai[n]*sizeof(int));
	for(i = 0; i < n; i ++)
	{
		m->aj[m->ai[i]] = i;
		signal[i] = 0;
		for(j = m->ai[i]+1; j < m->ai[i+1]; j ++)
		{
			int jdx = i + rand()%(n-i);
			while(!signal[jdx])
				jdx = i + rand()%(n-i);
			m->aj[j] = jdx;
			signal[jdx] = 0;
		}
		qsort(m->aj+m->ai[i], m->ai[i+1]-m->ai[i], sizeof(int), _cmp);
		for(j = m->ai[i]; j < m->ai[i+1]; j ++)
			signal[m->aj[j]] = 1;
	}
	//Init a
	m->a  = malloc(m->ai[n]*sizeof(double));
	for(i = 0; i < n; i ++)
	{
		m->a[m->ai[i]] = 1.0;
		for(j = m->ai[i]+1; j < m->ai[i+1]; j ++)
		{
			m->a[j] = (double)(rand())/(double)(RAND_MAX);
		}
	}
	//Init nz
	m->nz = signal;
	for(i = 0; i < n; i ++) m->nz[i] = 0;
	for(i = 0; i < n; i ++)
		for(j = m->ai[i]; j < m->ai[i+1]; j ++)
			m->nz[m->aj[j]] ++;
	return 0;
}

int Init_Diagonal(Matrix *m, int n)
{
	int i;
	m->n = n;
	m->ai = malloc((n+1)*sizeof(int));
	m->aj = malloc(n*sizeof(int));
	m->nz = malloc(n*sizeof(int));
	m->a  = malloc(n*sizeof(double));
	m->ai[0] = 0;
	for(i = 0; i < n; i ++)
	{
		m->ai[i+1] = m->ai[i] + 1;
		m->aj[i]   = i;
		m->nz[i]   = 1;
		m->a[i]    = 1.1;
	}
}

int Init_3D7(Matrix *m, int M, int N, int L, int max)
{
	int i,j,k;
	int n = m->n = M*N*L;
	int *aaj = NULL;
	m->ai = malloc((n+1)*sizeof(int));
	//Init ai
	m->ai[0] = 0;
	for(i = 0; i < M; i ++)
	{
		for(j = 0; j < N; j ++)
		{
			for(k = 0; k < L; k ++)
			{
				int nz = 1;
				if(i < M-1) nz ++;
				if(j < N-1) nz ++;
				if(k < L-1) nz ++;
				m->ai[i*N*L+j*L+k+1] = m->ai[i*N*L+j*L+k] + nz;
			}
		}
	}
	//Init aj
	m->aj = malloc(m->ai[n]*sizeof(int));
	aaj = m->aj;
	for(i = 0; i < M; i ++)
	{
		for(j = 0; j < N; j ++)
		{
			for(k = 0; k < L; k ++)
			{
				*aaj++ = i*N*L+j*L+k;
				if(k < L-1) *aaj++ = i*N*L+j*L+k+1; 
				if(j < N-1) *aaj++ = i*N*L+j*L+k+L; 
				if(i < M-1) *aaj++ = i*N*L+j*L+k+N*L;
			}
		}
	}
	//Init a
	m->a  = malloc(m->ai[n]*sizeof(double));
	for(i = 0; i < n; i ++)
	{
		m->a[m->ai[i]] = 0.01;
		for(j = m->ai[i]+1; j < m->ai[i+1]; j ++)
		{
			m->a[j] = (double)(rand())/(double)(RAND_MAX);
		}
	}
	//Init nz
	m->nz = malloc(n*sizeof(int));;
	for(i = 0; i < n; i ++) m->nz[i] = 0;
	for(i = 0; i < n; i ++)
		for(j = m->ai[i]; j < m->ai[i+1]; j ++)
			m->nz[m->aj[j]] ++;
	return 0;
}

int Init_3D13(Matrix *m, int M, int N, int L, int max)
{
	int i,j,k;
	int n = m->n = M*N*L;
	int *aaj = NULL;
	m->ai = malloc((n+1)*sizeof(int));
	//Init ai
	m->ai[0] = 0;
	for(i = 0; i < M; i ++)
	{
		for(j = 0; j < N; j ++)
		{
			for(k = 0; k < L; k ++)
			{
				int nz = 1;
				if(i < M-1) nz ++;
				if(i < M-2) nz ++;
				if(j < N-1) nz ++;
				if(j < N-2) nz ++;
				if(k < L-1) nz ++;
				if(k < L-2) nz ++;
				m->ai[i*N*L+j*L+k+1] = m->ai[i*N*L+j*L+k] + nz;
			}
		}
	}
	//Init aj
	m->aj = malloc(m->ai[n]*sizeof(int));
	aaj = m->aj;
	for(i = 0; i < M; i ++)
	{
		for(j = 0; j < N; j ++)
		{
			for(k = 0; k < L; k ++)
			{
				*aaj++ = i*N*L+j*L+k;
				if(k < L-1) *aaj++ = i*N*L+j*L+k+1; 
				if(k < L-2) *aaj++ = i*N*L+j*L+k+2; 
				if(j < N-1) *aaj++ = i*N*L+j*L+k+L; 
				if(j < N-2) *aaj++ = i*N*L+j*L+k+2*L; 
				if(i < M-1) *aaj++ = i*N*L+j*L+k+N*L;
				if(i < M-2) *aaj++ = i*N*L+j*L+k+2*N*L;
			}
		}
	}
	//Init a
	m->a  = malloc(m->ai[n]*sizeof(double));
	for(i = 0; i < n; i ++)
	{
		m->a[m->ai[i]] = 0.01;
		for(j = m->ai[i]+1; j < m->ai[i+1]; j ++)
		{
			m->a[j] = (double)(rand())/(double)(RAND_MAX);
		}
	}
	//Init nz
	m->nz = malloc(n*sizeof(int));;
	for(i = 0; i < n; i ++) m->nz[i] = 0;
	for(i = 0; i < n; i ++)
		for(j = m->ai[i]; j < m->ai[i+1]; j ++)
			m->nz[m->aj[j]] ++;
	return 0;
}

int Output_Matrix(Matrix *m)
{
	int n = m->n;
	int A[50][50];
	if(n > 32)
	{
		printf("Matrix is too large %d\n", n);
		return -1;
	}
	int i,j,k;
	printf("   ");
	for(i = 0; i < n; i ++)
	{
		printf("%2d ", m->ai[i+1]-m->ai[i]);
	}
	printf("\n====================================================================================================\n");
	for(i = 0; i < n; i ++)
		for(j = 0; j < n; j ++)
			A[i][j] = 0;
	for(i = 0; i < n; i ++)
		for(j = m->ai[i]; j < m->ai[i+1]; j ++)
			A[m->aj[j]][i] = 1;
	for(i = 0; i < n; i ++)
	{
		printf("%2d|", m->nz[i]);
		for(j = 0; j < n; j ++)
		{
			if(A[i][j] == 1)
				printf(" 1 ");
			else
				printf("   ");
		}
		printf("\n");
	}
	printf("====================================================================================================\n");
	for(i = 0; i < n; i ++)
	{
		for(j = m->ai[i]; j < m->ai[i+1]; j ++)
		{
			printf("%2d ", m->aj[j]);
		}
		printf("\n");
	}
	printf("====================================================================================================\n");
	return 0;
}

int Finalize_Matrix(Matrix *m)
{
	if(m->a  != NULL) free(m->a);
	if(m->ai != NULL) free(m->ai);
	if(m->aj != NULL) free(m->aj);
	if(m->nz != NULL) free(m->nz);
	return 0;
}

int Finalize_Matrix_Reorder(Matrix_Reorder *m)
{
	if(m->a   != NULL) free(m->a);
	if(m->ai  != NULL) free(m->ai);
	if(m->aj  != NULL) free(m->aj);
	if(m->m2r != NULL) free(m->m2r);
	if(m->r2m != NULL) free(m->r2m);
	return 0;
}

int Finalize_Matrix_Special(Matrix_Special *s)
{
	int i;
	if(s->nz     != NULL) free(s->nz);
	if(s->ai     != NULL) free(s->ai);
	if(s->aj     != NULL) free(s->aj);
	if(s->ll     != NULL) free(s->ll);
	if(s->ths    != NULL) free(s->ths);
	if(s->ll_ths != NULL) free(s->ll_ths);
	if(s->aa     != NULL) free(s->aa);
	for(i = 0; i < 32; i ++)
	{
		if(s->cnt[i] != NULL) free(s->cnt[i]);
	}

}
int Finalize_Matrix_WXL(Matrix_WXL *m)
{
	int i;
	if(m->ai    != NULL) free(m->ai);
	if(m->aj    != NULL) free(m->aj);
	if(m->aa    != NULL) free(m->aa);
	if(m->mai    != NULL) free(m->ai);
	if(m->maj    != NULL) free(m->aj);
	if(m->maa    != NULL) free(m->aa);
	if(m->idx   != NULL) free(m->idx);
	if(m->level != NULL) free(m->level);
	for(i = 0; i < 32; i ++)
		if(m->cnt[i] != NULL) free(m->cnt[i]);
}
int Finalize_Matrix_Compress(Matrix_Compress *m)
{
	int i;
	if(m->ai    != NULL) free(m->ai);
	if(m->aj    != NULL) free(m->aj);
	if(m->aa    != NULL) free(m->aa);
	if(m->mai    != NULL) free(m->ai);
	if(m->maj    != NULL) free(m->aj);
	if(m->maa    != NULL) free(m->aa);
	if(m->_idx  != NULL) free(m->_idx);
	if(m->level != NULL) free(m->level);
	for(i = 0; i < 32; i ++)
		if(m->cnt[i] != NULL) free(m->cnt[i]);
}
int Finalize_Matrix_Compress_General(Matrix_Compress_General *m)
{
	int i;
	if(m->ai    != NULL) free(m->ai);
	if(m->aj    != NULL) free(m->aj);
	if(m->aa    != NULL) free(m->aa);
	if(m->mai    != NULL) free(m->ai);
	if(m->maj    != NULL) free(m->aj);
	if(m->maa    != NULL) free(m->aa);
	if(m->_idx  != NULL) free(m->_idx);
	if(m->level != NULL) free(m->level);
	if(m->plevel != NULL) free(m->plevel);
	for(i = 0; i < CONSUMERS; i ++)
		if(m->cnt[i] != NULL) free(m->cnt[i]);
}

int Init_Random_Vector(double *x, int N)
{
	int i;
	for(i = 0; i < N; i ++) x[i] = 100.0;//(double)(rand())/(double)(RAND_MAX);
}

int Init_Zero_Vector(double *x, int N)
{
	int i;
	for(i = 0; i < N; i ++) x[i] = 0;
}

int Basic_Solver_Compress(double *x, Matrix_Compress *M, double *b)
{
	int l, id, i, j, k;
	if(x != b) for(i = 0; i < M->N; i ++) x[i] = b[i];
	int levels = M->levels;
	int n      = M->n;
	int begin  = 0;
	for(l = 0; l < levels; l ++)
	{
		for(k = 0; k < 9; k ++)
		{
			for(id = 0; id < 32; id ++)
			{
				for(j = M->idx[id][9*l+k]; j < M->idx[id][9*l+k+1]; j ++)
				{
					x[M->aj[j]] -= M->aa[j] * x[M->ai[j]];
				}
			}
		}
	}
	return 0;
}

int Basic_Solver_Compress_General(double *x, Matrix_Compress_General *M, double *b)
{
	int l, id, i, j, k;
	if(x != b) for(i = 0; i < M->N; i ++) x[i] = b[i];
	int levels = M->plevels;
	int n      = M->n;
	int begin  = 0;
	for(l = 0; l < levels; l ++)
	{
		for(k = 0; k < 9; k ++)
		{
			for(id = 0; id < 32; id ++)
			{
				for(j = M->idx[id][9*l+k]; j < M->idx[id][9*l+k+1]; j ++)
				{
					x[M->aj[j]] -= M->aa[j] * x[M->ai[j]];
				}
			}
		}
	}
	return 0;
}

int Basic_Solver_WXL(double *x, Matrix_WXL *M, double *b)
{
	int l,i,j,k,id;
	int N   = M->N;
	int nnz = M->nnz;
	if(x != b) for(i = 0; i < N; i ++) x[i] = b[i];
	//for(i = 0; i < nnz; i ++)
	//{
	//	x[M->aj[i]] -= M->aa[i] * x[M->ai[i]];
	//}
	int begin = 0;
	for(l = 0; l < M->levels; l ++)
	{
		int nn  = M->level[l];
		int dnn = nn & 0x1F;
		int pnn = nn >> 5;
		int offset = 32*(begin+nn);
		int *idx = M->idx + begin*9;
		for(k = 0; k < 9; k ++)
		{
			for(id = 0; id < 32; id ++)
			{
				int len  = pnn + (id < dnn ? 1 : 0);
				int left = id*pnn + (id < dnn ? id : dnn);
				for(i = 0; i < len; i ++)
				{
					for(j = idx[k*nn+left+i]; j < idx[k*nn+left+i+1]; j ++)
					{
						x[M->aj[j]] -= M->aa[j] * x[M->ai[j]];
					}
				}
			}
		}
		begin += nn;
	}
}

int Basic_Solver_CSC(double *x, double *a, int *ai, int *aj, double *b, int N)
{
	int i,j;
	if(x != b) for(i = 0; i < N; i ++) x[i] = b[i];
	for(i = 0; i < N; i ++)
	{
		//if(i == 0) printf("%8d %.20lf %.20lf\n", i, a[ai[i]], x[i]);
        x[i] *= a[ai[i]];
		for(j = ai[i]+1; j < ai[i+1]; j ++)
			x[aj[j]] -= a[j] * x[i];
	}
	return 0;
}
int Analyse_Matrix(Matrix *m)
{
	int i,j,k;
	int n = m->n;
	int *level = malloc(n*sizeof(int));
	int *count = malloc(n*sizeof(int));
	int levels = 0; 
	for(i = 0; i < n; i ++) level[i] = 0;
	for(i = 0; i < n; i ++) count[i] = 0;
	for(i = 0; i < n; i ++)
	{
		for(j = m->ai[i]+1; j < m->ai[i+1]; j ++)
		{
			int ll = level[i] + 1;
			if(ll > level[m->aj[j]]) level[m->aj[j]] = ll;
		}
	}
	for(i = 0; i < n; i ++)
	{
		count[level[i]] ++;
		if(level[i] > levels) levels = level[i];
	}
	levels ++;
	for(i = 0; i < levels; i ++)
	{
		//printf("%4d: %6d\n", i, count[i]);
	}
	free(level);
	free(count);
	return 0;
}

int Load_From_File(Matrix *m, char *name)
{
	FILE *file = fopen(name, "rb");
	int i,j,k;
	int N, nnz, exp;
	fread(&N,   sizeof(int), 1, file);
	fread(&nnz, sizeof(int), 1, file);
    printf("File Name: %s, #rows=#cols=%d, #nonzeros=%d\n", name, N, nnz);
	if((N%VECTOR) != 0)
    {
        exp = VECTOR - (N%VECTOR);
    }
    m->n = N+exp;
	m->ai = malloc((N+1+exp)*sizeof(int));
	fread(m->ai, sizeof(int), N+1, file);
	m->aj = malloc((m->ai[N]+exp)*sizeof(int));
	m->a  = malloc((m->ai[N]+exp)*sizeof(double));
	fread(m->aj, sizeof(int),    m->ai[N], file);
	fread(m->a , sizeof(double), m->ai[N], file);
	m->nz = NULL;
    if(exp > 0)
    {
        for(i = 0; i < exp; i ++)
            m->ai[N+1+i] = m->ai[N+i] + 1;
        for(i = 0; i < exp; i ++)
            m->aj[m->ai[N+i]] = N+i;
        for(i = 0; i < exp; i ++)
            m->a [m->ai[N+i]] = 1.0;
    }
	//for(i = 0; i < N; i ++)
	//{
	//	if(i != m->aj[m->ai[i]]) printf("Error 1\n");
	//	for(j = m->ai[i]+1; j < m->ai[i+1]; j ++)
	//	{
	//		if(m->aj[j] < i)
	//			printf("Error 2\n");
	//	}
	//}
	fclose(file);
	return N+exp;
}

void Output_Result(char *_name, int SF, int SB, int MF, int MB, double R)
{
	char name[200];
	strcpy(name, _name);
	char *n = name + strlen(name);
	while(*n != '.')
	{
		*n = 0;
		n --;
	}
	*n = 0;
	n --;
	while(*n != '/')
		n --;
	n ++;

	double ratio = ((double)SF)/((double)MF);

	printf("GET RESULT: SLAVE(MFlops %8d MB/s %8d) MASTER(MFlops %8d MB/s %8d) Ratio: %.6lf R/T: %9.3lf% Name %s\n", SF, SB, MF, MB, ratio, R*100, n);
	//printf("%s %d %d\n", n, F, B);
	return;
}

int main(int argc, char **argv)
{
	athread_init();
    printf("======================================\n");	
	int STEPS = 100, s;
	int N = 10;
	int mnz = 5;
	double *b = NULL;
	Matrix m;
	Matrix_Reorder  r;
	//Matrix_Compress wxl;
	Matrix_Compress_General wxl;

	N = Load_From_File(&m, argv[1]);
	long long size = (sizeof(double)+sizeof(int))*m.ai[N] + sizeof(double)*2*N + sizeof(int)*N;
	int flop = m.ai[N] * 2 - N;
	double parallelism = Reorder(&r,  &m);
    printf("The parallelism is %.6lf\n", parallelism);
    printf("PRODUCER_CONSUMER_ROWS:%d CACHE_X_V4:%d\n", PRODUCER_CONSUMER_ROWS, CACHE_X_V4);
	Finalize_Matrix(&m);
	Compress_General(&wxl, &r);

	double *x0 = malloc(N*sizeof(double));
	double *x1 = malloc(N*sizeof(double));
	if(b == NULL)
	{
		b  = malloc(N*sizeof(double));
		Init_Random_Vector(b, N);
	}
	else
	{
		Init_Random_Vector(b, N);
	}
	
	struct timeval t1, t2;
    printf("======Test SwSpTRSV %d times=======\n", STEPS);
    gettimeofday(&t1, NULL);
    Parallel_Solver_Compress_General(x1, &wxl, b, STEPS);
    gettimeofday(&t2, NULL);
    printf("Filename: %s, PRODUCER_CONSUMER_ROWS: %d CACHE_X_V4: %d Average time is %.6lfs, Average MFlops is %.6lf\n", argv[1], PRODUCER_CONSUMER_ROWS, CACHE_X_V4, TIME(t1,t2)/STEPS, flop*1e-6/(TIME(t1,t2)/STEPS));
	
	Basic_Solver_CSC(x0, r.a, r.ai, r.aj, b, N);
	
	printf("=========Check================================\n");
	//check(x0, x1, N);
	special_check(x0, x1, N);
	printf("=========Check Over===========================\n");

	free(x0);
	free(x1);
	free(b);
	Finalize_Matrix_Reorder(&r);
	Finalize_Matrix_Compress_General(&wxl);
    printf("======================================\n");	
	return 0;
}
