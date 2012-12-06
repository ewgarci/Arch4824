/*
 * Matrices are represented with struct smat.
 * Do not roll your own struct for representing matrices! Use smat.
 *
 * You are welcome to use any of the helper functions exported
 * below to allocate/deallocate smat's.
 * smat_printf can be very useful for debugging (along with gdb).
 */
#ifndef _MATRIX_HELPERS_H_
#define _MATRIX_HELPERS_H_

#include <stddef.h>

/*
 * Our representation of a matrix.
 * We call it "Sommelier" Matrix (smat for short)
 *
 * Note that vectors are a particular case in which cols == 1,
 * and scalars have rows == cols == 1.
 */
struct smat {
	double **data;
	unsigned int rows;
	unsigned int cols;
};

struct smat *smat_calloc(size_t rows, size_t cols);
void smat_free(struct smat *mat);
void smat_printf(const char *title, const struct smat *mat);

/*
 * You probably won't need the functions below; you're welcome to use them
 * as long as in doing that you're not duplicating smat's functionality.
 */
void free_2d_double(double **matrix);
void *calloc_2d_double(size_t rows, size_t cols);


#endif /* _MATRIX_HELPERS_H_ */
