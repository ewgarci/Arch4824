/*
 * Optimized Unicore implementations of the matrix algorithms.
 * Note that the result matrix r is not initisialised to 0.
 */
#include <stdio.h>
#include <string.h>
#include "strassen.c"
#include "smat.h"

void smat_memset(struct smat *m)
{
	int i;
	for (i = 0; i < m->rows; i++)
		memset (m->data[i],'0', m->cols * sizeof(double));
}

/* Matrix multiplication, r = a x b */
void smat_mult(const struct smat *a, const struct smat *b, struct smat *r)
{
	
	smat_memset(r);

	strassen(a->data, b->data, r->data, a->rows);

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
