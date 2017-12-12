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
    R->level = max_level;
    R->levels = malloc(R->level*sizeof(int));
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
    for(i = 0; i < max_level; i ++) R->levels[i] = cnt[i+1]-cnt[i];
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

	//Check Reorder
#if 0
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

void CSC2CSR(Matrix_Reorder *csr, Matrix_Reorder *csc)
{
    int level = csr->level = csc->level;
    int n     = csr->n     = csc->n;
    int nnz   = csc->ai[n]-csc->ai[0];
    int m = n;
    csr->ai = malloc((m+1)*sizeof(int));
    csr->aj = malloc(nnz*sizeof(int));
    csr->a  = malloc(nnz*sizeof(double));
    int i, j;
    for(i = 0; i < m+1; i ++) csr->ai[i] = 0;
    for(i = 0; i < m; i ++)
        for(j = csc->ai[i]; j < csc->ai[i+1]; j ++) csr->ai[csc->aj[j]+1] ++;
    for(i = 0; i < m; i ++) csr->ai[i+1] += csr->ai[i];
    int    **ajp = malloc(m*sizeof(int*));
    double **aap = malloc(m*sizeof(double*));
    for(i = 0; i < m; i ++)
    {
        ajp[i] = csr->aj + csr->ai[i];
        aap[i] = csr->a  + csr->ai[i];
    }
    for(i = 0; i < m; i ++)
        for(j = csc->ai[i]; j < csc->ai[i+1]; j ++)
        {
            *ajp[csc->aj[j]]++ = i;
            *aap[csc->aj[j]]++ = csc->a [j];
        }
    for(i = 0; i < m; i ++)
    {
        if(ajp[i] != csr->aj + csr->ai[i+1]) printf("Error in %d\n", i);
        if(aap[i] != csr->a  + csr->ai[i+1]) printf("Error in %d\n", i);
    }

    csr->levels = malloc(level*sizeof(int));
    csr->m2r = malloc(n*sizeof(int));
    csr->r2m = malloc(n*sizeof(int));

    memcpy(csr->levels, csc->levels, level*sizeof(int));
    memcpy(csr->m2r,    csc->m2r,    n*sizeof(int));
    memcpy(csr->r2m,    csc->r2m,    n*sizeof(int));
    free(ajp);
    free(aap);
    return;
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
    if(m->levels != NULL)   free(m->levels);
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
	if(m->_cnt != NULL) free(m->_cnt);
}

int Init_Vector_100(double *x, int N)
{
	int i;
	for(i = 0; i < N; i ++) x[i] = 100.0;//(double)(rand())/(double)(RAND_MAX);
}

int Init_Zero_Vector(double *x, int N)
{
	int i;
	for(i = 0; i < N; i ++) x[i] = 0;
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
    printf("======================================\n");	
	int (*Compress_General)(Matrix_Compress_General*, Matrix_Reorder*) = NULL;
    int (*Parallel_Solver_Compress_General)(double*, Matrix_Compress_General*, double*, int) = NULL; 
    athread_init();
	int STEPS = 100, s;
	int N = 10;
	int mnz = 5;
	double *b = NULL;
	Matrix m;
	Matrix_Reorder  r, r2;
	Matrix_Compress_General wxl;
    
	N = Load_From_File(&m, argv[1]);
	long long size = (sizeof(double)+sizeof(int))*m.ai[N] + sizeof(double)*2*N + sizeof(int)*N;
	int flop = m.ai[N] * 2 - N;
	double parallelism = Reorder(&r,  &m);
	Finalize_Matrix(&m);


    printf("The parallelism is %.6lf\n", parallelism);
    if(parallelism < 32)
    {
        Compress_General = Compress_General_R2B64;
        Parallel_Solver_Compress_General = Parallel_Solver_Compress_General_R2B64;
    }
    else if(parallelism < 1024)
    {
        Compress_General = Compress_General_R4B512;
        Parallel_Solver_Compress_General = Parallel_Solver_Compress_General_R4B512;
    }
    else
    {
        Compress_General = Compress_General_R8B1024;
        Parallel_Solver_Compress_General = Parallel_Solver_Compress_General_R8B1024;
    }

	Compress_General(&wxl, &r);

	double *x0 = malloc(N*sizeof(double));
	double *x1 = malloc(N*sizeof(double));
	if(b == NULL)
	{
		b  = malloc(N*sizeof(double));
		Init_Vector_100(b, N);
	}
	else
	{
		Init_Vector_100(b, N);
	}
	
	struct timeval t1, t2;

    printf("======Test SwSpTRSV %d times=======\n", STEPS);
    gettimeofday(&t1, NULL);
    Parallel_Solver_Compress_General(x1, &wxl, b, STEPS);
    gettimeofday(&t2, NULL);
    printf("SwSpTRSV: Filename is %s Average time is %.6lfs, Average MFlops is %.6lf\n", argv[1], TIME(t1,t2)/STEPS, flop*1e-6/(TIME(t1,t2)/STEPS));
	Finalize_Matrix_Compress_General(&wxl);
    
    printf("======Test Serial %d times =======\n", STEPS);
    gettimeofday(&t1, NULL);
    for(s = 0; s < STEPS; s ++)
    {
        Basic_Solver_CSC(x0, r.a, r.ai, r.aj, b, N);
    }
    gettimeofday(&t2, NULL);
    printf("Serial: Filename is %s Average time is %.6lfs, Average MFlops is %.6lf\n", argv[1], TIME(t1,t2)/STEPS, flop*1e-6/(TIME(t1,t2)/STEPS));
    
    printf("======Test level-sets %d times =======\n", STEPS);
    CSC2CSR(&r2, &r);
    gettimeofday(&t1, NULL);
    Parallel_Solver_CSR(x1, r2.a, r2.ai, r2.aj, r2.levels, b, N, r2.level, STEPS); 
    gettimeofday(&t2, NULL);
    printf("level-sets: Filename is %s Average time is %.6lfs, Average MFlops is %.6lf\n", argv[1], TIME(t1,t2)/STEPS, flop*1e-6/(TIME(t1,t2)/STEPS));
	Finalize_Matrix_Reorder(&r2);
	Finalize_Matrix_Reorder(&r);
    printf("======================================\n");	

	free(x0);
	free(x1);
	free(b);
	return 0;
}
