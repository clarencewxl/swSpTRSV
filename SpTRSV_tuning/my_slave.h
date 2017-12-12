/*************************************************************************
	> File Name: my_slave.h
	> Author: Wang Xinliang
	> mail: clarencewxl@gmail.com
	> Created Time: Tue 11 Oct 2016 03:34:45 PM CST
 ************************************************************************/

#ifndef _MY_SLAVE_H_
#define _MY_SLAVE_H_

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

#define SLAVES 64
#define SLAVESLOG2 6
#define SLAVES2NUM 0x3F
#define COLS 8
#define COLSLOG2 3
#define COLS2NUM 0x07
#define ROWS 8
#define ROWSLOG2 3
#define ROWS2NUM 0x07
#define COL(x) (x & 0x07)
#define ROW(x) ((x & 0x38) >> 3)
#define REG_PUTR(var, dst) asm volatile ("putr %0,%1\n"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1\n"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0\n":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0\n":"=r"(var))
#define ROWSYN  athread_syn(ROW_SCOPE,0xff)
#define COLSYN  athread_syn(COL_SCOPE,0xff)
#define ALLSYN  athread_syn(ARRAY_SCOPE,0xffff)
#define WXL_DMA_SET(d,mode,size,reply) {\
	  dma_set_size(d, size); \
	  dma_set_bsize(d, size); \
	  dma_set_stepsize(d, 0); \
	  dma_set_op(d, mode); \
	  dma_set_mode(d, PE_MODE); \
	  dma_set_reply(d, reply); \
}
#define WXL_DMA_SET_NOSIZE(d,mode,reply) {\
	  dma_set_stepsize(d, 0); \
	  dma_set_op(d, mode); \
	  dma_set_mode(d, PE_MODE); \
	  dma_set_reply(d, reply); \
}
#define WXL_DMA_SET_SIZE(d,size) {\
	  dma_set_size(d, (size)); \
	  dma_set_bsize(d, (size)); \
}
#define WXL_DMA(d,l,r,size) {WXL_DMA_SET_SIZE(&d,size); dma(d,(long)(l),(long)(r));}
#define WXL_DMA_NEW(d,l,r,size,count) {WXL_DMA_SET_SIZE(&d,size); dma(d,(long)(l),(long)(r)); count ++;}

typedef struct
{
	long long *data;
	long long *buf;
	int M,N,L;
} Info;

typedef struct
{
	char *dst, *src;
	int size;
} Memcpy_Info;

typedef struct
{
	double *x, *a, *b;
	int *ai, *aj, *nz, *levels, *finish_count;
	int size;
} Solver_Info;

#endif
