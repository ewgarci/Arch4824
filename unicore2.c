/*
 * Optimized Unicore implementations of the matrix algorithms. Traverses Matrix Multiplication in columns
 * Note that the result matrix r is not initisialised to 0.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "smat.h"

void *calloc_2d_double(size_t rows, size_t cols);

void smat_memset(const struct smat *m)
{
	int i;
	for (i = 0; i < m->rows; i++)
		memset (m->data[i],'0', m->cols * sizeof(double));
}


/* Matrix multiplication, r = a x b */ 
void smat_mult(const struct smat *a, const struct smat *b, struct smat *r)
{
	int i, j, k;
	
	smat_memset(r);
	/*
	for (i = 0; i < r->rows; i++) {
		for (k = 0; k < r->rows; k++) {
			for (j = 0; j < r->rows; j++) {
				r->data[i][j] += a->data[i][k] * b->data[k][j];
			}
		}
	}
	
	double tmp[r->rows][r->rows];
	for (i = 0; i < r->rows; ++i)
		for (j = 0; j < r->rows; ++j)
			tmp[i][j] =  b->data[j][i];
		for (i = 0; i < r->rows; ++i)
			for (j = 0; j < r->rows; ++j)
				for (k = 0; k < r->rows; ++k)
					r->data[i][j] += a->data[i][k] * tmp[j][k];
	*/
	
		define SM (32 / sizeof (double))
		for (i = 0; i < r->rows; i += SM)
			for (j = 0; j < r->rows; j += SM)
				for (k = 0; k < r->rows; k += SM)
					for (i2 = 0, rres = &r->data[i][j],	rmul1 = &a->data[i][k]; i2 < SM; ++i2, rres += r->rows, rmul1 += r->rows)
							for (k2 = 0, rmul2 = &b->data[k][j]; k2 < SM; ++k2, rmul2 += r->rows)
		for (j2 = 0; j2 < SM; ++j2)
					rres[j2] += rmul1[k2] * rmul2[j2];
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

void *calloc_2d_double(size_t rows, size_t cols)
{
	double *flat;
	double **mem;
	int i;

	flat = calloc(rows, cols * sizeof(double));
	if (flat == NULL) {
		fprintf(stderr, "%s: Error allocating %Zu bytes: ",
			__func__, rows * cols * sizeof(double));
		perror(NULL);
		exit(EXIT_FAILURE);
	}

	mem = malloc(rows * sizeof(double *));
	if (mem == NULL) {
		perror("calloc_2d_double");
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < rows; i++)
		mem[i] = flat + i * cols;

	return mem;
}
