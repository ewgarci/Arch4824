/*
 * Array allocation/deallocation helpers.
 * NOTE: matrices are allocated row-by-row in contiguous virtual address space.
 */
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

#include "smat.h"

void free_2d_double(double **matrix)
{
	if (matrix == NULL)
		return;
	if (matrix[0])
		free(matrix[0]);
	free(matrix);
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
/*
void smat_memset(struct smat *m)
{
	int i;
	
	for (i = 0; i < m->rows; i++)
		memset (m->data[i],'0', m->cols * sizeof(double));
		
	
}
*/
/* Print array in an octave-friendly format */
void smat_printf(const char *title, const struct smat *mat)
{
	int is_scalar = mat->rows == 1 && mat->cols == 1;
	int i, j;

	if (title)
		printf("# name: %s\n", title);
	if (is_scalar) {
		printf("# type: scalar\n");
	} else {
		printf("# type: matrix\n");
		printf("# rows: %u\n", mat->rows);
		printf("# columns: %u\n", mat->cols);
	}
	for (i = 0; i < mat->rows; i++) {
		for (j = 0; j < mat->cols; j++)
			printf(" %.17e", mat->data[i][j]);
		printf("\n");
	}
}

void smat_free(struct smat *mat)
{
	if (mat == NULL)
		return;

	free_2d_double(mat->data);
	free(mat);
}

struct smat *smat_calloc(size_t rows, size_t cols)
{
	struct smat *m;

	m = malloc(sizeof(struct smat));
	if (m == NULL) {
		perror("calloc_smat");
		exit(EXIT_FAILURE);
	}

	m->rows = rows;
	m->cols = cols;
	m->data = calloc_2d_double(rows, cols);

	return m;
}

