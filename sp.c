/*
 * 
 * Note that the result matrix r is not initisialised to 0.
 */
#include <stdio.h>
#include <string.h>
#include "smat.h"
#include "mm.c"
#include "strassenParallel.c"



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
	strassenParallel(a->data, b->data, r->data, a->rows);

}

/*
 * Matrix-vector multiplication, r = a x v
 * NOTE: a is the matrix, v is the vector. This means v->cols is always 1.
 * The result is always a vector (r->cols == 1).
 */
void smat_vect(const struct smat *a, const struct smat *v, struct smat *r)
{
	parallelOperation(a->data, v->data, r->data, a->rows, VECTOR);
}

/* Matrix addition, i.e. r = a + b */
void smat_add(const struct smat *a, const struct smat *b, struct smat *r)
{
	parallelOperation(a->data, b->data, r->data, a->rows, ADD);
}

/* Scale matrix a by constant alpha */
void smat_scale(const struct smat *a, const struct smat *alpha, struct smat *r)
{
	parallelOperation(a->data, alpha->data, r->data, a->rows, SCALE);
}
