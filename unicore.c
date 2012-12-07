/*
 * Optimized Unicore implementations of the matrix algorithms. Uses Strassen Method for Multiplication 
 * Note that the result matrix r is not initisialised to 0.
 */
#include <stdio.h>
#include <string.h>
#include "strassen.c"
#include "smat.h"

void *calloc_2d_double(size_t rows, size_t cols);

/**********************************************************************************************
* Reproduced from
* http://stackoverflow.com/questions/9908917/optimized-static-padding-for-strassens-odd-matrices
**********************************************************************************************/

/*	
 *	Here Q is equal to the BREAK value in the strassen
 *	multiplication. 
 */
int get_best_pud_up_value(int actual_size, int Q) {
    int cnt = 0;
    int n = actual_size;
    while(n > Q) {
        cnt++;
        n /= 2;
    }

	/*
    result should be smallest value such that:
    result >= actual_size AND
    result % (1<<cnt) == 0
	*/
	
    if (actual_size % (1<<cnt) == 0) {
        return actual_size;
    } else {
        return actual_size + (1<<cnt) - actual_size % (1<<cnt);
    }
}

/*	In order for Strassen's Method to work for all matricies
 *	we find the next power of 2 sized matrix and pad the extra 
 *	rows and columns with 0's. As an optimization to the 
 *	code we can get away with padding to the next m * 2^k matrix.
 */
void *calloc_2d_double(size_t rows, size_t cols)
{
	double *flat;
	double **mem;
	int i;
	int N = get_best_pud_up_value(rows, 16);

	flat = calloc(N, N * sizeof(double));
	if (flat == NULL) {
		fprintf(stderr, "%s: Error allocating %Zu bytes: ",
			__func__, N * N * sizeof(double));
		perror(NULL);
		exit(EXIT_FAILURE);
	}

	mem = malloc(N * sizeof(double *));
	if (mem == NULL) {
		perror("calloc_2d_double");
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < N; i++)
		mem[i] = flat + i * N;

	return mem;
}


/* Matrix multiplication, r = a x b */
void smat_mult(const struct smat *a, const struct smat *b, struct smat *r)
{

	strassen(a->data, b->data, r->data, get_best_pud_up_value(r->rows, 16));
//	printf("strassen num = %d\n", r->strassenN);

}

/*****************************************************************
 * Reproduced from http://stackoverflow.com/questions/12289235/simple-and-fast-matrix-vector-multiplication-in-c-c
 ******************************************************************/
double vectors_dot_prod(double *x, double **y, int n)
{
    double res = 0.0;
    int i = 0;
    for (; i <= n-4; i+=4)
    {
		res += (x[i] * y[i][0] +
                x[i+1] * y[i+1][0] +
                x[i+2] * y[i+2][0] +
                x[i+3] * y[i+3][0]);
    }
    for (; i < n; i++)
    {
        res += x[i] * y[i][0];
    }
    return res;
}

/*
 * Matrix-vector multiplication, r = a x v
 * NOTE: a is the matrix, v is the vector. This means v->cols is always 1.
 * The result is always a vector (r->cols == 1).
 */
void smat_vect(const struct smat *a, const struct smat *v, struct smat *r)
{
	int i;

    for (i = 0; i < a->rows; i++)
       r->data[i][0] = vectors_dot_prod(a->data[i], v->data, a->cols);
		
}
	
/* Matrix addition, i.e. r = a + b */
void smat_add(const struct smat *a, const struct smat *b, struct smat *r)
{
	int i, j;

	for (i = 0; i < a->rows; i++)
		for (j = 0; j < a->cols; j++)
			r->data[i][j] = a->data[i][j] + b->data[i][j];
}

/* Scale matrix a by constant alpha */
void smat_scale(const struct smat *a, const struct smat *alpha, struct smat *r)
{
	double factor = alpha->data[0][0];
	int i, j;

	for (i = 0; i < a->rows; i++)
		for (j = 0; j < a->cols; j++)
			r->data[i][j] = a->data[i][j] * factor;
}



