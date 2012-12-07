#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include "job_api.h"

#define CPUCORES	8		
#define MULTIPLY	0
#define ADD			1
#define SCALE		2
#define VECTOR		3


double **aa, **bb, **cc;
int matrixSize;
void parallelOperation(double **a, double **b, double **c, int size, int operation);

static void mm(void *args)
{
	struct work_struct *s = args;
	int i,j,k;

	//printf("[%d]MM size = %d, threads = %d \n", s->id, matrixSize, CPUCORES);

	for (i = s->id; i < matrixSize; i += CPUCORES) {
		for (k = 0; k < matrixSize; k++) {
				for (j = 0; j < matrixSize; j++) {
					cc[i][j] += aa[i][k] * bb[k][j];
					//printf("\t[%d] cc[%d][%d]= %lf, aa= %lf, bb= %lf\n", k, i, j, cc[i][j], aa[i][k], bb[k][j]);	
			}
			//printf("cc[%d][%d]= %lf\n", i, j, cc[i][j]);
		}
	}
	//printf("%d", s->id);
	if (s->id)
		job_exit();
  
}

/*****************************************************************
 * Reproduced from http://stackoverflow.com/questions/12289235/simple-and-fast-matrix-vector-multiplication-in-c-c
 ******************************************************************/

static double vectors_dot_prod(double *x, double **y)
{
    double res = 0.0;
    int i = 0;
	//printf("[1]res= %lf, x[%d]= %lf, y[%d][0]= %lf\n", res, i, x[i], i, y[i][0]);
    for (; i <= matrixSize-4; i+=4)
    {
		//printf("[2]res= %lf, x[%d]= %lf, y[%d][0]= %lf\n", res, i, x[i], i, y[i][0]);
		res += (x[i] * y[i][0] +
                x[i+1] * y[i+1][0] +
                x[i+2] * y[i+2][0] +
                x[i+3] * y[i+3][0]);
    }
    for (; i < matrixSize; i++)
    {
        res += x[i] * y[i][0];
    }
    return res;
	
}

static void mv(void *args)
{
	struct work_struct *s = args;
	int i;

	//printf("[%d]MV size = %d, threads = %d \n", s->id, matrixSize, CPUCORES);

	for (i = s->id; i < matrixSize; i += CPUCORES){
		cc[i][0] = vectors_dot_prod(aa[i], bb);
	}
	//printf("%d", s->id);
	 if (s->id)
		job_exit();
}

static void add(void *args)
{
	struct work_struct *s = args;
	int i,j;

	//printf("[%d]ADD size = %d, threads = %d \n", s->id, matrixSize, CPUCORES);

	for (i = s->id; i < matrixSize; i += CPUCORES) {
		for (j = 0; j < matrixSize; j++) {
			cc[i][j] = aa[i][j] + bb[i][j];
		}
	}  
	//printf("%d", s->id);
	if(s->id)
		job_exit();
}

static void scale(void *args)
{
	struct work_struct *s = args;
	double factor = bb[0][0];
	int i,j;

	//printf("[%d]SCALE size = %d, threads = %d \n", s->id, matrixSize, CPUCORES);

	for (i = s->id; i < matrixSize; i += CPUCORES) {
		for (j = 0; j < matrixSize; j++) {
			cc[i][j] = aa[i][j] * factor;
		}
	}
	//printf("%d", s->id);
	if (s->id)
		job_exit();
}

 void parallelOperation(double **a, double **b, double **c, int size, int operation)
 {
	int i;
	struct work_struct *jobs[CPUCORES];

	aa = a;
	bb = b;
	cc = c;

	matrixSize = size;

	for (i = 0; i < CPUCORES; i++) {
		jobs[i] = malloc(sizeof(struct work_struct));
		if (jobs[i] == NULL) {
			perror(NULL);
			exit(1);
		}
		jobs[i]->id = i;
	}
	
	job_init();
	
	switch(operation){
		case MULTIPLY:{
	//		printf("Multiply[%d]: ", size);
			for (i = 1; i < CPUCORES; i++)
				job_create(mm, jobs[i], 0);
			mm(jobs[0]);
		}break;

		case VECTOR:{
	//		printf("Vector[%d]: ", size);
			for (i = 1; i < CPUCORES; i++)
				job_create(mv, jobs[i], 0);
			mv(jobs[0]);
		}break;

		case ADD:{
	//		printf("Add[%d]: ", size);
			for (i = 1; i < CPUCORES; i++)
				job_create(add, jobs[i], 0);
			add(jobs[0]);
		}break;

		case SCALE:{
	//		printf("Scale[%d]: ", size);
			for (i = 1; i < CPUCORES; i++)
				job_create(scale, jobs[i], 0);
			scale(jobs[0]);
		}break;
	}
		

	for (i = 1; i < CPUCORES; i++)
		job_join(jobs[i]);

//	printf("\n");
}

