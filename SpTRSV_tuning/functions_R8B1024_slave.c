/*************************************************************************
  > File Name: functions_slave.c
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

static inline int EMPTY_GETR(void) { int var; asm volatile ("rcsr %0,0x34\n":"=r"(var)); return var; }
static inline int EMPTY_GETC(void) { int var; asm volatile ("rcsr %0,0x35\n":"=r"(var)); return var; }
static inline int FULL_PUT(void) { int var; asm volatile ("rcsr %0,0x36\n":"=r"(var)); return var; }

static inline uint64_t slave_get_rtc()
{
	unsigned long rpcc;
	asm volatile ("rcsr %0, 4":"=&r"(rpcc)::"memory");
	return rpcc;
}

#define REG_GETR_I2(v0,v1,var)     REG_GETR(var); v0=((long long*)(&(var)))[0]; v1=((long long*)(&(var)))[1]
#define REG_GETC_I2(v0,v1,var)     REG_GETC(var); v0=((long long*)(&(var)))[0]; v1=((long long*)(&(var)))[1]
#define REG_PUTR_I2(v0,v1,var,dst) ((long long*)(&(var)))[0]=(v0); ((long long*)(&(var)))[1]=(v1); REG_PUTR(var,dst)
#define REG_PUTC_I2(v0,v1,var,dst) ((long long*)(&(var)))[0]=(v0); ((long long*)(&(var)))[1]=(v1); REG_PUTC(var,dst)
#define REG_GETR_DI(v0,v1,var)	   REG_GETR(var); v0=((double*)(&(var)))[0]; v1=((int*)(&(var)))[2]
#define REG_GETC_DI(v0,v1,var)     REG_GETC(var); v0=((double*)(&(var)))[0]; v1=((int*)(&(var)))[2]
#define REG_PUTR_DI(v0,v1,var,dst) ((double*)(&(var)))[0]=(v0); ((int*)(&(var)))[2]=(v1); REG_PUTR(var,dst)
#define REG_PUTC_DI(v0,v1,var,dst) ((double*)(&(var)))[0]=(v0); ((int*)(&(var)))[2]=(v1); REG_PUTC(var,dst)

void wxl_dma_wait(int *v, int C)
{
	dma_wait(v,C);
	return;
}

static void Consumer(doublev4 *x, doublev4 *b, int *_level, int *_cnt, int levels, int n)
{
	int rid = ROW(_MYID);
	int cid = COL(_MYID);
	int  id = rid * CONSUMER_COLS + cid;
	int tgt = cid + CONSUMER_COLS;
	int  tn = (n>>CONSUMERS_LOG2) + (id < (n&CONSUMERS_2NUM) ? 1 : 0);
	volatile int reply = 0;
	volatile int COUNT = 0;
	dma_desc dma_get;
	dma_desc dma_put;
	doublev4 buffer[CONSUMER_BUFFER_V4];
	doublev4 _value[CACHE_X_V4];
	double *value = _value;
	int   cnt[CONSUMER_CACHE];
	int level[CONSUMER_CACHE];
	
    WXL_DMA_SET_NOSIZE(&dma_get, DMA_GET, &reply);
	WXL_DMA_SET_NOSIZE(&dma_put, DMA_PUT, &reply);
	
	int begin_global = 0;
	int begin_local  = 0;
	int dlen         = 0;
	int buf_idx      = CONSUMER_BUFFER_V4;
	//Init value
	int nn           = _level[0];
	int len          = (nn>>CONSUMERS_LOG2) + (id < (nn&CONSUMERS_2NUM) ? 1 : 0);
	begin_global    += nn;
	begin_local     += len;
	tn              -= len;
	b               += len;
	dlen             = MIN(tn, CACHE_X_V4 - (begin_local&CACHE_X_V4_2NUM));
	if(dlen > 0)
	{
		WXL_DMA_NEW(dma_get, b, &_value[len&CACHE_X_V4_2NUM], dlen*sizeof(doublev4), COUNT);
		tn    -= dlen;
		b     += dlen;

	}
	dlen = MIN(tn, begin_local&CACHE_X_V4_2NUM);
	if(dlen > 0)
	{
		WXL_DMA_NEW(dma_get, b, &_value[0], dlen*sizeof(doublev4), COUNT);
		tn    -= dlen;
		b     += dlen;

	}
	wxl_dma_wait(&reply, COUNT);

	if(buf_idx == CONSUMER_BUFFER_V4 && tn > 0)
	{
		buf_idx = 0;
		dlen    = MIN(tn, CONSUMER_BUFFER_V4);
		WXL_DMA_NEW(dma_get, b, buffer, dlen*sizeof(doublev4), COUNT);
		tn    -= dlen;
		b     += dlen;

	}
	
	int i,ii,j,k,iid;
	for(i = 0, ii = 0; i < levels-1; i ++, ii = (ii+1)&CONSUMER_CACHE_2NUM)
	{
		if(ii == 0)
		{
			//Init level and cnt
			int tmp = MIN(CONSUMER_CACHE, levels-i);
			WXL_DMA_NEW(dma_get, _level+i+1, level, tmp*sizeof(int), COUNT);
			WXL_DMA_NEW(dma_get, _cnt+i,     cnt,   tmp*sizeof(int), COUNT);
			wxl_dma_wait(&reply, COUNT);

		}
		nn  = level[ii];
		iid = (id + CONSUMERS - (begin_global&CONSUMERS_2NUM)) & CONSUMERS_2NUM;
		len = (nn>>CONSUMERS_LOG2) + (iid < (nn&CONSUMERS_2NUM) ? 1 : 0);
		//Get data from Producer
		int  cntt = cnt[ii];
		int _cntt = cntt - (cntt&VECTOR_2NUM);
		for(j = 0; j < _cntt; j += VECTOR)
		{
			doublev4 vvv0, vvv1, vvv2, vvv3;
			REG_GETR(vvv0);
			REG_GETR(vvv1);
			REG_GETR(vvv2);
			REG_GETR(vvv3);
			int aj0 = ((int*)(&vvv0))[2];
			int aj1 = ((int*)(&vvv1))[2];
			int aj2 = ((int*)(&vvv2))[2];
			int aj3 = ((int*)(&vvv3))[2];
			value[((aj0>>CONSUMERS_LOG2)&(CACHE_X_V4_2NUM<<VECTOR_LOG2))+(aj0&VECTOR_2NUM)] -= ((double*)(&vvv0))[0];
			value[((aj1>>CONSUMERS_LOG2)&(CACHE_X_V4_2NUM<<VECTOR_LOG2))+(aj1&VECTOR_2NUM)] -= ((double*)(&vvv1))[0];
			value[((aj2>>CONSUMERS_LOG2)&(CACHE_X_V4_2NUM<<VECTOR_LOG2))+(aj2&VECTOR_2NUM)] -= ((double*)(&vvv2))[0];
			value[((aj3>>CONSUMERS_LOG2)&(CACHE_X_V4_2NUM<<VECTOR_LOG2))+(aj3&VECTOR_2NUM)] -= ((double*)(&vvv3))[0];
		}
		for(j; j < cntt; j ++)
		{
			doublev4 vvv;
			REG_GETR(vvv);
			int aj = ((int*)(&vvv))[2];
			value[((aj>>CONSUMERS_LOG2)&(CACHE_X_V4_2NUM<<VECTOR_LOG2))+(aj&VECTOR_2NUM)] -= ((double*)(&vvv))[0];
		}	
		//Put data to Producer
		dlen         = MIN(len, CACHE_X_V4 - (begin_local&CACHE_X_V4_2NUM));
		int _dlen    =   dlen - (dlen&VECTOR_2NUM);
		doublev4 *vv = _value + (begin_local&CACHE_X_V4_2NUM);
		for(j = 0; j < _dlen; j += VECTOR)
		{
			REG_PUTR(vv[j+0], tgt);
			REG_PUTR(vv[j+1], tgt);
			REG_PUTR(vv[j+2], tgt);
			REG_PUTR(vv[j+3], tgt);
		}
		for(j; j < dlen; j ++)
		{
			REG_PUTR(vv[j+0], tgt);
		}
		int tmp = len-dlen;
		dlen    = MAX(0, tmp);
		_dlen   = dlen - (dlen&VECTOR_2NUM); 
		for(j = 0; j < _dlen; j += VECTOR)
		{
			REG_PUTR(_value[j+0], tgt);
			REG_PUTR(_value[j+1], tgt);
			REG_PUTR(_value[j+2], tgt);
			REG_PUTR(_value[j+3], tgt);
		}
		for(j; j < dlen; j ++)
		{
			REG_PUTR(_value[j+0], tgt);
		}

		//Supply data from buffer
		dlen = MIN(len, CACHE_X_V4 - (begin_local&CACHE_X_V4_2NUM));
		vv   = _value + (begin_local&CACHE_X_V4_2NUM);
		while(dlen > 0 && buf_idx < CONSUMER_BUFFER_V4)
		{
			dma_wait(&reply, COUNT);
			int ddlen = MIN(dlen, (CONSUMER_BUFFER_V4 - buf_idx));
			for(j = 0; j < ddlen; j ++)
				*vv++ = buffer[buf_idx+j];
			buf_idx += ddlen;
			dlen    -= ddlen;
			if(buf_idx == CONSUMER_BUFFER_V4 && tn > 0)
			{
				buf_idx = 0;
				int dddlen = MIN(tn, CONSUMER_BUFFER_V4);
				WXL_DMA_NEW(dma_get, b, buffer, dddlen*sizeof(doublev4), COUNT);
				tn     -= dddlen;
				b      += dddlen;

			}
		}
		dlen = len - MIN(len, CACHE_X_V4 - (begin_local&CACHE_X_V4_2NUM));
		vv   = _value;
		while(dlen > 0 && buf_idx < CONSUMER_BUFFER_V4)
		{
			dma_wait(&reply, COUNT);
			int ddlen = MIN(dlen, (CONSUMER_BUFFER_V4 - buf_idx));
			for(j = 0; j < ddlen; j ++)
				*vv++ = buffer[buf_idx+j];
			buf_idx += ddlen;
			dlen    -= ddlen;
			if(buf_idx == CONSUMER_BUFFER_V4 && tn > 0)
			{
				buf_idx = 0;
				int dddlen = MIN(tn, CONSUMER_BUFFER_V4);
				WXL_DMA_NEW(dma_get, b, buffer, dddlen*sizeof(doublev4), COUNT);
				tn     -= dddlen;
				b      += dddlen;

			}
		}
		begin_global += nn;
		begin_local  += len;
	}	
	return;
}

static void Producer_General(doublev4 *x, doublev4 *b, double *_aa, unsigned int *_ai, unsigned int *_aj, int *_idx, int *_level, int levels, int n)
{
	int rid = ROW(_MYID);
	int cid = COL(_MYID) - CONSUMER_COLS;
	int  id = rid * PRODUCER_COLS + cid;
	int tgt = (rid+1)&PRODUCER_ROWS_2NUM;
	volatile int reply = 0;
	volatile int COUNT = 0;
	dma_desc dma_get;
	dma_desc dma_put;
	WXL_DMA_SET_NOSIZE(&dma_get, DMA_GET, &reply);
	WXL_DMA_SET_NOSIZE(&dma_put, DMA_PUT, &reply);
	
	doublev4 *x_original = x;
	
    //Value
	doublev4   _value[CACHE_X_V4];
	double     *value                                                 = _value;
	doublev4 (*_value2)[CACHE_X_V4>>PRODUCER_ROWS_LOG2]               = _value;
	double    (*value2)[CACHE_X_V4>>(PRODUCER_ROWS_LOG2-VECTOR_LOG2)] = _value;

	//Buffer for mass memory access
	doublev4 buffer[PRODUCER_BUFFER_V4];

	//Matrix A
	double       aa[PRODUCER_BUFFER_A];
	unsigned int aj[PRODUCER_BUFFER_A];
    
	//Index & Idx
	int    level[PRODUCER_CACHE_LEVEL];
	int    idx  [PRODUCER_CACHE_IDX*MAX_SUB_LEVELS];

	int idx00        = _idx[0];
	int idxnn        = _idx[levels*(PRODUCER_ROWS+1)];
	int an           =  idxnn-idx00;
	int buf_idx      = 0;
	int   a_idx      = 0;
	_aa             += idx00;
	_aj             += idx00;
	int len          = MIN(PRODUCER_BUFFER_A, an);
	unsigned int begin_global = 0;

    if(len > 0)
	{
		WXL_DMA_NEW(dma_get, _aa, aa, len*sizeof(double),       COUNT);
		WXL_DMA_NEW(dma_get, _aj, aj, len*sizeof(unsigned int), COUNT);
		_aa += len;
		_aj += len;
		an  -= len;

	}
	int i,ii,iii,j,k;
	int last_idx = idx00;
	for(i = 0, ii = 0, iii = 0; i < levels; i ++, ii = (ii+1)&PRODUCER_CACHE_IDX_2NUM, iii = (iii+1)&PRODUCER_CACHE_LEVEL_2NUM)
	{
		uint64_t tt0, tt1;
		tt0 = slave_get_rtc();

		if(ii == 0)
		{
			int idx_len = MIN(PRODUCER_CACHE_IDX, levels-i);
			WXL_DMA_NEW(dma_get, _idx + MAX_SUB_LEVELS*i + 1, idx, idx_len*MAX_SUB_LEVELS*sizeof(int), COUNT);
			wxl_dma_wait(&reply, COUNT);

		}
		if(iii == 0)
		{
			int idx_len = MIN(PRODUCER_CACHE_LEVEL, levels-i);
			WXL_DMA_NEW(dma_get, _level+i, level, idx_len*sizeof(int), COUNT);
			wxl_dma_wait(&reply, COUNT);

		}
		int nn   = level[iii];
		if(nn <= 0)
		{
			nn = -nn;
			int len  = CACHE_X_V4;
			unsigned int base = begin_global<<VECTOR_LOG2; 
			WXL_DMA_NEW(dma_get, x_original+CACHE_X_V4*nn, _value, len*sizeof(doublev4), COUNT);

			dma_wait(&reply, COUNT);
			last_idx = idx[ii*MAX_SUB_LEVELS];
			for(k = 1; k < MAX_SUB_LEVELS; k ++)
			{
				//exchange data
				if(k > 1)
				{
					for(j = 0; j < CACHE_X_V4; j += VECTOR)
					{
						REG_PUTC(_value[j+0], tgt);
						REG_PUTC(_value[j+1], tgt);
						REG_PUTC(_value[j+2], tgt);
						REG_PUTC(_value[j+3], tgt);
						REG_GETC(_value[j+0]);
						REG_GETC(_value[j+1]);
						REG_GETC(_value[j+2]);
						REG_GETC(_value[j+3]);
					}
				}
				//compute
				int didx = idx[ii*MAX_SUB_LEVELS + k] - last_idx;
				while(didx > 0)
				{
					dma_wait(&reply, COUNT);
					int dlen = MIN(didx, PRODUCER_BUFFER_A - a_idx);
					int _dlen = dlen - (dlen&VECTOR_2NUM);
					for(j = 0; j < _dlen; j += VECTOR)
					{
						doublev4 vvv0, vvv1, vvv2, vvv3;
						unsigned int dst0 = (((aj[a_idx+j+0]&BIG_BLOCK_2NUM))>>VECTOR_LOG2)&CONSUMER_COLS_2NUM;
						unsigned int dst1 = (((aj[a_idx+j+1]&BIG_BLOCK_2NUM))>>VECTOR_LOG2)&CONSUMER_COLS_2NUM;
						unsigned int dst2 = (((aj[a_idx+j+2]&BIG_BLOCK_2NUM))>>VECTOR_LOG2)&CONSUMER_COLS_2NUM;
						unsigned int dst3 = (((aj[a_idx+j+3]&BIG_BLOCK_2NUM))>>VECTOR_LOG2)&CONSUMER_COLS_2NUM;
						((unsigned int*)(&vvv0))[2] = (aj[a_idx+j+0]&BIG_BLOCK_2NUM)+base;
						((unsigned int*)(&vvv1))[2] = (aj[a_idx+j+1]&BIG_BLOCK_2NUM)+base;
						((unsigned int*)(&vvv2))[2] = (aj[a_idx+j+2]&BIG_BLOCK_2NUM)+base;
						((unsigned int*)(&vvv3))[2] = (aj[a_idx+j+3]&BIG_BLOCK_2NUM)+base;
						((double*)(&vvv0))[0] = value[(aj[a_idx+j+0]>>BIG_BLOCK_LOG2)&CACHE_X_2NUM] * aa[a_idx+j+0];
						((double*)(&vvv1))[0] = value[(aj[a_idx+j+1]>>BIG_BLOCK_LOG2)&CACHE_X_2NUM] * aa[a_idx+j+1];
						((double*)(&vvv2))[0] = value[(aj[a_idx+j+2]>>BIG_BLOCK_LOG2)&CACHE_X_2NUM] * aa[a_idx+j+2];
						((double*)(&vvv3))[0] = value[(aj[a_idx+j+3]>>BIG_BLOCK_LOG2)&CACHE_X_2NUM] * aa[a_idx+j+3];
						REG_PUTR(vvv0, dst0);
						REG_PUTR(vvv1, dst1);
						REG_PUTR(vvv2, dst2);
						REG_PUTR(vvv3, dst3);
					}
					for(j; j < dlen; j ++)
					{
						doublev4 vvv0;
						unsigned int dst0 = (((aj[a_idx+j+0]&BIG_BLOCK_2NUM))>>VECTOR_LOG2)&CONSUMER_COLS_2NUM;
						((unsigned int*)(&vvv0))[2] = (aj[a_idx+j+0]&BIG_BLOCK_2NUM)+base;
						((double*)(&vvv0))[0] = value[(aj[a_idx+j+0]>>BIG_BLOCK_LOG2)&CACHE_X_2NUM] * aa[a_idx+j+0];
						REG_PUTR(vvv0, dst0);
					}
					a_idx += dlen;
					didx -= dlen;
					if(a_idx == PRODUCER_BUFFER_A && an > 0)
					{
						a_idx   = 0;
						int dan = MIN(PRODUCER_BUFFER_A, an);
						WXL_DMA_NEW(dma_get, _aa, aa, dan*sizeof(double), COUNT);
						WXL_DMA_NEW(dma_get, _aj, aj, dan*sizeof(unsigned int),    COUNT);
						_aa += dan;
						_aj += dan;
						an  -= dan;

					}
				}
				last_idx = idx[ii*MAX_SUB_LEVELS+k];
			}
		}
		else
		{
			unsigned int base = (begin_global+nn)<<VECTOR_LOG2; 
			int iid  = (id + PRODUCERS -(begin_global&PRODUCERS_2NUM))&PRODUCERS_2NUM;
			int len  = (nn>>PRODUCERS_LOG2) + (iid < (nn&PRODUCERS_2NUM) ? 1 : 0);
			//Level 0
			if(len > 0)
			{
				if(i == 0)
				{
					WXL_DMA_NEW(dma_get, b, _value, len*sizeof(doublev4), COUNT);
					dma_wait(&reply, COUNT);

				}
				else
				{
					int _len = len - (len&VECTOR_2NUM);
					for(j = 0; j < _len; j += VECTOR)
					{
						REG_GETR(_value[j+0]);
						REG_GETR(_value[j+1]);
						REG_GETR(_value[j+2]);
						REG_GETR(_value[j+3]);
					}
					for(j; j < len; j ++)
					{
						REG_GETR(_value[j+0]);
					}
				}
				//compute
				int didx = idx[ii*MAX_SUB_LEVELS]-last_idx;
				while(didx > 0)
				{
					dma_wait(&reply, COUNT);
					int dlen = MIN(didx, PRODUCER_BUFFER_A - a_idx);
					for(j = 0; j < dlen; j ++)
					{
						unsigned int aai =  aj[a_idx+j] >> BIG_BLOCK_LOG2;
                        unsigned int aaj =  aj[a_idx+j]  & BIG_BLOCK_2NUM;
                        value[aaj] = (aai == aaj ? (value[aaj] * aa[a_idx+j]) : (value[aaj] - value[aai] * aa[a_idx+j]));
                        //value[aaj] = value[aaj] - value[aai] * aa[a_idx+j];
					}
					a_idx += dlen;
					didx -= dlen;
					if(a_idx == PRODUCER_BUFFER_A && an > 0)
					{
						a_idx = 0;
						int dan = MIN(PRODUCER_BUFFER_A, an);
						WXL_DMA_NEW(dma_get, _aa, aa, dan*sizeof(double), COUNT);
						WXL_DMA_NEW(dma_get, _aj, aj, dan*sizeof(unsigned int),    COUNT);
						_aa += dan;
						_aj += dan;
						an  -= dan;

					}
				}
				last_idx = idx[ii*MAX_SUB_LEVELS];
				//Write Data to main memory
				//if((CONSUMER_BUFFER_V4 - buf_idx) >= len)
				if((CONSUMER_BUFFER_V4 - buf_idx) >= len && ((x - x_original + buf_idx + len)&CACHE_X_V4_2NUM) > ((x - x_original)&CACHE_X_V4_2NUM))
				{
					for(j = 0; j < len; j ++)
					{
						buffer[buf_idx+j] = _value[j];
					}
					buf_idx += len;
				}
				else
				{
					if(buf_idx > 0)
					{
						WXL_DMA_NEW(dma_put, x, buffer, buf_idx*sizeof(doublev4), COUNT);
						x      += buf_idx;

						buf_idx = 0;
					}
					WXL_DMA_NEW(dma_put, x, _value, len*sizeof(doublev4), COUNT);
					x     += len;

					dma_wait(&reply, COUNT);
				}
			}
			//Level 1-9
			if(nn <= CACHE_X)
			{
				int max_dlen = (nn + PRODUCERS_2NUM)>>PRODUCERS_LOG2;
				int _max_dlen = max_dlen - (max_dlen&VECTOR_2NUM);
				for(k = 0; k < PRODUCER_ROWS - 1; k ++)
				{
					for(j = 0; j < _max_dlen; j += VECTOR)
					{
						REG_PUTC(_value2[k][j+0], tgt);
						REG_PUTC(_value2[k][j+1], tgt);
						REG_PUTC(_value2[k][j+2], tgt);
						REG_PUTC(_value2[k][j+3], tgt);
						REG_GETC(_value2[k+1][j+0]);
						REG_GETC(_value2[k+1][j+1]);
						REG_GETC(_value2[k+1][j+2]);
						REG_GETC(_value2[k+1][j+3]);
					}	
					for(j; j < max_dlen; j ++)
					{
						REG_PUTC(_value2[k][j+0], tgt);
						REG_GETC(_value2[k+1][j+0]);
					}
				}
                int didx = idx[ii*MAX_SUB_LEVELS + PRODUCER_ROWS] - last_idx;
                while(didx > 0)
                {
                    dma_wait(&reply, COUNT);
                    int dlen = MIN(didx, PRODUCER_BUFFER_A-a_idx);
                    int _dlen = dlen - (dlen&VECTOR_2NUM);
                    for(j = 0; j < _dlen; j += VECTOR)
                    {
                        doublev4 vvv0, vvv1, vvv2, vvv3;
                        unsigned int dst0 = (((aj[a_idx+j+0]&BIG_BLOCK_2NUM)+base)>>VECTOR_LOG2)&CONSUMER_COLS_2NUM;
                        unsigned int dst1 = (((aj[a_idx+j+1]&BIG_BLOCK_2NUM)+base)>>VECTOR_LOG2)&CONSUMER_COLS_2NUM;
                        unsigned int dst2 = (((aj[a_idx+j+2]&BIG_BLOCK_2NUM)+base)>>VECTOR_LOG2)&CONSUMER_COLS_2NUM;
                        unsigned int dst3 = (((aj[a_idx+j+3]&BIG_BLOCK_2NUM)+base)>>VECTOR_LOG2)&CONSUMER_COLS_2NUM;
                        ((unsigned int*)(&vvv0))[2] = (aj[a_idx+j+0]&BIG_BLOCK_2NUM)+base;
                        ((unsigned int*)(&vvv1))[2] = (aj[a_idx+j+1]&BIG_BLOCK_2NUM)+base;
                        ((unsigned int*)(&vvv2))[2] = (aj[a_idx+j+2]&BIG_BLOCK_2NUM)+base;
                        ((unsigned int*)(&vvv3))[2] = (aj[a_idx+j+3]&BIG_BLOCK_2NUM)+base;
                        ((double*)(&vvv0))[0] = value2[aj[a_idx+j+0]>>(BIG_BLOCK_LOG2+CACHE_X_LOG2)][(aj[a_idx+j+0]>>BIG_BLOCK_LOG2)&CACHE_X_2NUM] * aa[a_idx+j+0];
                        ((double*)(&vvv1))[0] = value2[aj[a_idx+j+1]>>(BIG_BLOCK_LOG2+CACHE_X_LOG2)][(aj[a_idx+j+1]>>BIG_BLOCK_LOG2)&CACHE_X_2NUM] * aa[a_idx+j+1];
                        ((double*)(&vvv2))[0] = value2[aj[a_idx+j+2]>>(BIG_BLOCK_LOG2+CACHE_X_LOG2)][(aj[a_idx+j+2]>>BIG_BLOCK_LOG2)&CACHE_X_2NUM] * aa[a_idx+j+2];
                        ((double*)(&vvv3))[0] = value2[aj[a_idx+j+3]>>(BIG_BLOCK_LOG2+CACHE_X_LOG2)][(aj[a_idx+j+3]>>BIG_BLOCK_LOG2)&CACHE_X_2NUM] * aa[a_idx+j+3];
                        REG_PUTR(vvv0, dst0);
                        REG_PUTR(vvv1, dst1);
                        REG_PUTR(vvv2, dst2);
                        REG_PUTR(vvv3, dst3);
                    }
                    for(j; j < dlen; j ++)
                    {
                        doublev4 vvv0;
                        unsigned int dst0 = (((aj[a_idx+j+0]&BIG_BLOCK_2NUM)+base)>>VECTOR_LOG2)&CONSUMER_COLS_2NUM;
                        ((unsigned int*)(&vvv0))[2] = (aj[a_idx+j+0]&BIG_BLOCK_2NUM)+base;
                        ((double*)(&vvv0))[0] = value2[aj[a_idx+j+0]>>(BIG_BLOCK_LOG2+CACHE_X_LOG2)][(aj[a_idx+j+0]>>BIG_BLOCK_LOG2)&CACHE_X_2NUM] * aa[a_idx+j+0];
                        REG_PUTR(vvv0, dst0);
                    }
                    a_idx += dlen;
                    didx -= dlen;
                    if(a_idx == PRODUCER_BUFFER_A && an > 0)
                    {
                        a_idx = 0;
                        int dan = MIN(PRODUCER_BUFFER_A, an);
                        WXL_DMA_NEW(dma_get, _aa, aa, dan*sizeof(double), COUNT);
                        WXL_DMA_NEW(dma_get, _aj, aj, dan*sizeof(unsigned int),    COUNT);
                        _aa += dan;
                        _aj += dan;
                        an  -= dan;

                    }
                }
                last_idx = idx[ii*MAX_SUB_LEVELS + PRODUCER_ROWS];
			}
			else
			{
				for(k = 1; k < MAX_SUB_LEVELS; k ++)
				{
					//exchange data
					if(k > 1)
					{
                        int max_dlen = (nn + PRODUCERS_2NUM)>>PRODUCERS_LOG2;
                        int _max_dlen = max_dlen - (max_dlen&VECTOR_2NUM);
						for(j = 0; j < _max_dlen; j += VECTOR)
						{
							REG_PUTC(_value[j+0], tgt);
							REG_PUTC(_value[j+1], tgt);
							REG_PUTC(_value[j+2], tgt);
							REG_PUTC(_value[j+3], tgt);
							REG_GETC(_value[j+0]);
							REG_GETC(_value[j+1]);
							REG_GETC(_value[j+2]);
							REG_GETC(_value[j+3]);
						}
						for(j; j < max_dlen; j ++)
						{
							REG_PUTC(_value[j], tgt);
							REG_GETC(_value[j]);
						}
					}
					//compute
					int didx = idx[ii*MAX_SUB_LEVELS+k]-last_idx;
					while(didx > 0)
					{
						dma_wait(&reply, COUNT);
						int dlen = MIN(didx, PRODUCER_BUFFER_A - a_idx);
						int _dlen = dlen - (dlen&VECTOR_2NUM);
						for(j = 0; j < _dlen; j += VECTOR)
						{
							doublev4 vvv0, vvv1, vvv2, vvv3;
							unsigned int dst0 = (((aj[a_idx+j+0]&BIG_BLOCK_2NUM)+base)>>VECTOR_LOG2)&CONSUMER_COLS_2NUM;
							unsigned int dst1 = (((aj[a_idx+j+1]&BIG_BLOCK_2NUM)+base)>>VECTOR_LOG2)&CONSUMER_COLS_2NUM;
							unsigned int dst2 = (((aj[a_idx+j+2]&BIG_BLOCK_2NUM)+base)>>VECTOR_LOG2)&CONSUMER_COLS_2NUM;
							unsigned int dst3 = (((aj[a_idx+j+3]&BIG_BLOCK_2NUM)+base)>>VECTOR_LOG2)&CONSUMER_COLS_2NUM;
							((unsigned int*)(&vvv0))[2] = (aj[a_idx+j+0]&BIG_BLOCK_2NUM)+base;
							((unsigned int*)(&vvv1))[2] = (aj[a_idx+j+1]&BIG_BLOCK_2NUM)+base;
							((unsigned int*)(&vvv2))[2] = (aj[a_idx+j+2]&BIG_BLOCK_2NUM)+base;
							((unsigned int*)(&vvv3))[2] = (aj[a_idx+j+3]&BIG_BLOCK_2NUM)+base;
                            ((double*)(&vvv0))[0] = value[(aj[a_idx+j+0]>>BIG_BLOCK_LOG2)&CACHE_X_2NUM] * aa[a_idx+j+0];
                            ((double*)(&vvv1))[0] = value[(aj[a_idx+j+1]>>BIG_BLOCK_LOG2)&CACHE_X_2NUM] * aa[a_idx+j+1];
                            ((double*)(&vvv2))[0] = value[(aj[a_idx+j+2]>>BIG_BLOCK_LOG2)&CACHE_X_2NUM] * aa[a_idx+j+2];
                            ((double*)(&vvv3))[0] = value[(aj[a_idx+j+3]>>BIG_BLOCK_LOG2)&CACHE_X_2NUM] * aa[a_idx+j+3];
                            //((double*)(&vvv0))[0] = value[(aj[a_idx+j+0]>>BIG_BLOCK_LOG2)] * aa[a_idx+j+0];
                            //((double*)(&vvv1))[0] = value[(aj[a_idx+j+1]>>BIG_BLOCK_LOG2)] * aa[a_idx+j+1];
                            //((double*)(&vvv2))[0] = value[(aj[a_idx+j+2]>>BIG_BLOCK_LOG2)] * aa[a_idx+j+2];
                            //((double*)(&vvv3))[0] = value[(aj[a_idx+j+3]>>BIG_BLOCK_LOG2)] * aa[a_idx+j+3];
							REG_PUTR(vvv0, dst0);
							REG_PUTR(vvv1, dst1);
							REG_PUTR(vvv2, dst2);
							REG_PUTR(vvv3, dst3);
						}
						for(j; j < dlen; j ++)
						{
							doublev4 vvv0;
							unsigned int dst0 = (((aj[a_idx+j+0]&BIG_BLOCK_2NUM)+base)>>VECTOR_LOG2)&CONSUMER_COLS_2NUM;
							((unsigned int*)(&vvv0))[2] = (aj[a_idx+j+0]&BIG_BLOCK_2NUM)+base;
                            ((double*)(&vvv0))[0] = value[(aj[a_idx+j+0]>>BIG_BLOCK_LOG2)&CACHE_X_2NUM] * aa[a_idx+j+0];
                            //((double*)(&vvv0))[0] = value[(aj[a_idx+j+0]>>BIG_BLOCK_LOG2)] * aa[a_idx+j+0];
							REG_PUTR(vvv0, dst0);
						}
						a_idx += dlen;
						didx  -= dlen;
						if(a_idx == PRODUCER_BUFFER_A && an > 0)
						{
							a_idx = 0;
							int dan = MIN(PRODUCER_BUFFER_A, an);
							WXL_DMA_NEW(dma_get, _aa, aa, dan*sizeof(double), COUNT);
							WXL_DMA_NEW(dma_get, _aj, aj, dan*sizeof(unsigned int),    COUNT);
							_aa += dan;
							_aj += dan;
							an  -= dan;

						}
					}
					last_idx = idx[ii*MAX_SUB_LEVELS + k];
				}
			}
			begin_global += nn;
		}
		//ROWSYN;
		athread_syn(ROW_SCOPE,0xF0);
		tt1 = slave_get_rtc();
	}
	if(buf_idx > 0)
	{
		WXL_DMA_NEW(dma_put, x, buffer, buf_idx*sizeof(doublev4), COUNT);
		x += buf_idx;

	}
	dma_wait(&reply, COUNT);
}

void General_Solver_CSC_Multiple_MPE_R8B1024(void *_ptr)
{
	Matrix_Compress_General m;
	volatile int reply = 0;
	volatile int COUNT = 0;
	dma_desc dma_get;
	WXL_DMA_SET_NOSIZE(&dma_get, DMA_GET, &reply);
	WXL_DMA_NEW(dma_get, _ptr, &m, sizeof(Matrix_Compress_General), COUNT);
	wxl_dma_wait(&reply, COUNT);
	int s = 0;
    for(s = 0; s < m.times; s ++)
    {
        if(ROW(_MYID) < PRODUCER_CONSUMER_ROWS)
        {
            if(COL(_MYID) < CONSUMER_COLS)
            {
                int      n = m.n;
                int     id = ROW(_MYID)*CONSUMER_COLS + COL(_MYID);
                int offset = id*(n>>CONSUMERS_LOG2) + (id < (n&CONSUMERS_2NUM) ? id : (n&CONSUMERS_2NUM));
                Consumer(m.x+VECTOR*offset, m.b+VECTOR*offset, m.level, m.cnt[id], m.levels, m.n);
            }
            else
            {
                int      n = m.n;
                int     id = ROW(_MYID)*PRODUCER_COLS + COL(_MYID) - CONSUMER_COLS;
                int offset = id*(n>>PRODUCERS_LOG2) + (id < (n&PRODUCERS_2NUM) ? id : (n&PRODUCERS_2NUM));
                Producer_General(m.x+VECTOR*offset, m.b+VECTOR*offset, m.aa, m.ai, m.aj, m.idx[id], m.plevel, m.plevels, m.n);
            }
        }
        ALLSYN;
    }
	return;
}

