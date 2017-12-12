/*************************************************************************
  > File Name: functions_levelsets_slave.c
  > Author: Wang Xinliang
  > mail: clarencewxl@gmail.com
  > Created Time: Mon 13 Feb 2017 06:35:18 PM CST
 ************************************************************************/

#define PRODUCER_CONSUMER_ROWS      8
#define PRODUCER_CONSUMER_ROWS_LOG2 3
#define CACHE_X_V4                  1024
#define CACHE_X_V4_LOG2             10

#include <stdio.h>
#include <stdlib.h>
#include <slave.h>
#include <simd.h>
#include <dma.h>
#include "my_slave.h"
#include "my_solver.h"
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

static void Parallel_Levelsets_better(double *x, double *aa, int *ai, int *aj, int *levels, double *b, int N, int L)
{
    double laa[4096];
    int  laj[4096];
    int  lai[1025];
    int   ll[1024];
    int myid = _MYID;
    int num_threads = COLS * PRODUCER_CONSUMER_ROWS;
    int offset = 0;
    int level = 0;
	volatile int reply = 0;
	volatile int COUNT = 0;
	dma_desc dma_get;
	WXL_DMA_SET_NOSIZE(&dma_get, DMA_GET, &reply);
    for(level = 0; level < L; level ++)
    {
        if((level&0x3FF) == 0)
        {
			WXL_DMA_NEW(dma_get, levels+level,     ll,   MIN(1024, L-level)*sizeof(int), COUNT);
			wxl_dma_wait(&reply, COUNT);
        }

        int tasks = ll[level&0x3FF];
        int my_tasks  = tasks/num_threads + (myid < (tasks%num_threads) ? 1 : 0);
        int my_offset = myid * (tasks/num_threads) + (myid < (tasks%num_threads) ? myid : (tasks%num_threads));
        int i,j,k,ii;
        for(k = 0; k < my_tasks; k += 1024)
        {
            int len = MIN(1024, my_tasks-k);
            WXL_DMA_NEW(dma_get, ai+k+offset+my_offset, lai, (len+1)*sizeof(int), COUNT);
			wxl_dma_wait(&reply, COUNT);
            i  = offset + my_offset + k;
            ii = 0;
            while(i < offset + my_offset + k + len)
            {
                int count = 0;
                int _ai = lai[ii];
                while((ii+count+1) <= len && (lai[ii+count+1]-_ai) <= 4096) count ++;
                if(count == 0)
                {
                    int sub_offset = _ai;
                    int sub_len = lai[ii+1]-lai[ii];
                    double xx = b[i];
                    while(sub_len > 4096)
                    {
                        WXL_DMA_NEW(dma_get, aa+sub_offset, laa, sizeof(double)*4096, COUNT);
                        WXL_DMA_NEW(dma_get, aj+sub_offset, laj, sizeof(int)*4096, COUNT);
                        dma_wait(&reply, COUNT);
                        for(j = 0; j < 4096; j ++)
                            xx -= laa[j] * x[laj[j]];
                        sub_len    -= 4096;
                        sub_offset += 4096;
                    }
                    WXL_DMA_NEW(dma_get, aa+sub_offset, laa, sizeof(double)*sub_len, COUNT);
                    WXL_DMA_NEW(dma_get, aj+sub_offset, laj, sizeof(int)*sub_len, COUNT);
                    dma_wait(&reply, COUNT);
                    for(j = 0; j < sub_len-1; j ++)
                        xx -= laa[j] * x[laj[j]];
                    x[i] = xx*laa[j];

                    i ++;
                    ii ++;
                }
                else
                {
                    int iii;
                    int size = lai[ii+count]-_ai;
                    WXL_DMA_NEW(dma_get, aa+_ai, laa, sizeof(double)*size, COUNT);
                    WXL_DMA_NEW(dma_get, aj+_ai, laj, sizeof(int)*size, COUNT);
			        dma_wait(&reply, COUNT);
                    for(iii = 0; iii < count; iii ++)
                    {
                        double xx = b[iii+i];
                        for(j = lai[iii+ii]; j < lai[iii+ii+1]-1; j ++)
                            xx -= laa[j-_ai] * x[laj[j-_ai]];
                        x[iii+i] = xx * laa[j-_ai];
                    }
                    i += count;
                    ii += count;
                }

            }
        }
        offset += tasks;
        ALLSYN;
    }
}

void Levelsets(void *_ptr)
{
    Matrix_LS m = *((Matrix_LS*)_ptr);
    int s;
    for(s = 0; s < m.times; s ++)
    {
        if(ROW(_MYID) < PRODUCER_CONSUMER_ROWS)
        {
            Parallel_Levelsets_better(m.x, m.aa, m.ai, m.aj, m.levels, m.b, m.N, m.L);
        }
        ALLSYN;
    }
}

