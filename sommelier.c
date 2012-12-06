/*
 * sommelier -- a matrix cruncher for fluid viscosity measurement.
 *
 * Pass the -h flag for usage information.
 */
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "sommelier.h"
#include "smat.h"

#define SMALL_M		20
#define SMALL_N		5
#define MEDIUM_M	71
#define MEDIUM_N	2
#define LARGE_M		201
#define LARGE_N		2

static const char progname[] = "sommelier";

static unsigned int matrix_size = 100;
static unsigned int random_seed = 10;
static unsigned int n_matrices = 10;
static struct smat *matrices;
static struct smat *vectors;
static struct smat *alphas;
static int debug;

static void __usage(FILE *fp)
{
	fprintf(fp, "Usage: %s [OPTIONS]\n", progname);
	fprintf(fp, "  Computes a set of array operations in a sequence.\n");
	fprintf(fp, "  These operations are: array multiplication, ");
	fprintf(fp, "array-vector multiplication, array addition, ");
	fprintf(fp, "and array scaling.\n");
	fprintf(fp, "\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -d: Print Debug information.\n");
	fprintf(fp, "  -h: Display this help message.\n");
	fprintf(fp, "  -m: Set matrix size. (default: %u)\n", matrix_size);
	fprintf(fp, "  -n: Set number of matrices. (default: %u)\n",
		n_matrices);
	fprintf(fp, "  -s: Set initial randomness seed. (default: %u)\n",
		random_seed);
	fprintf(fp, "  -t: test to be run. Choose from small, medium, large\n");
	fprintf(fp, "    These tests are just a shortcut for particular ");
	fprintf(fp, "(m,n) pairs.\n");
	fprintf(fp, "\n");
	fprintf(fp, "Notes:\n");
	fprintf(fp, "  - The input sequence is generated using the ");
	fprintf(fp, "pseudo-random generator random(3). ");
	fprintf(fp, "This sequence is guaranteed to be deterministic across ");
	fprintf(fp, "runs for a given random seed\n");
}

static inline void usage_err(void)
{
	__usage(stderr);
}

static inline void usage(void)
{
	__usage(stdout);
}

static void parse_args(int argc, char *argv[])
{
	int c;

	for (;;) {
		c = getopt(argc, argv, "dhm:n:s:t:");
		if (c < 0)
			break;
		switch (c) {
		case 'd':
			debug = 1;
			break;
		case 'h':
			usage();
			exit(EXIT_SUCCESS);
		case 'm':
			matrix_size = strtoul(optarg, NULL, 0);
			break;
		case 'n':
			n_matrices = strtoul(optarg, NULL, 0);
			break;
		case 's':
			random_seed = strtol(optarg, NULL, 0);
			break;
		case 't':
			if (!strncmp(optarg, "small", 6)) {
				matrix_size	= SMALL_M;
				n_matrices	= SMALL_N;
			} else if (!strncmp(optarg, "medium", 7)) {
				matrix_size	= MEDIUM_M;
				n_matrices	= MEDIUM_N;
			} else if (!strncmp(optarg, "large", 6)) {
				matrix_size	= LARGE_M;
				n_matrices	= LARGE_N;
			}
			break;
		default:
			usage_err();
			exit(EXIT_FAILURE);
		}
	}
}

static void matrix_fill_in(struct smat *arr)
{
	int i, j;

	srandom(++random_seed);

	for (i = 0; i < arr->rows; i++) {
		for (j = 0; j < arr->cols; j++) {
			long randval = random();

			arr->data[i][j] = random();
			arr->data[i][j] = arr->data[i][j] / RAND_MAX * 2.0;
			if (randval & 1)
				arr->data[i][j] *= -1;
		}
	}
}

static struct smat *__create(unsigned int n, size_t rows, size_t cols)
{
	struct smat *array;
	int i;

	array = calloc(n, sizeof(struct smat));
	if (array == NULL) {
		perror("__create");
		exit(EXIT_FAILURE);
	}

	for (i = 0; i < n; i++) {
		struct smat *m = &array[i];

		m->rows = rows;
		m->cols = cols;
		m->data = calloc_2d_double(rows, cols);

		matrix_fill_in(m);
	}

	return array;
}

static inline void create_matrices(void)
{
	matrices = __create(n_matrices + 1, matrix_size, matrix_size);
}

static inline void create_vectors(void)
{
	vectors = __create(n_matrices, matrix_size, 1);
}

static inline void create_alphas(void)
{
	alphas = __create(n_matrices, 1, 1);
}

static void debug_op(const char *op, const char *a_title, const struct smat *a,
		     const char *b_title, const struct smat *b,
		     const char *r_title, const struct smat *r)
{
	if (!debug)
		return;

	printf("# op: %s\n", op);
	smat_printf(a_title, a);
	smat_printf(b_title, b);
	smat_printf(r_title, r);
	printf("# ---\n");
}

static inline void
__smat_mult(const struct smat *a, const struct smat *b, struct smat *r)
{
	smat_mult(a, b, r);
	debug_op("mult", "a", a, "b", b, "r", r);
}

static inline void
__smat_scale(const struct smat *a, const struct smat *alpha, struct smat *r)
{
	smat_scale(a, alpha, r);
	debug_op("scale", "a", a, "alpha", alpha, "r", r);
}

static inline void
__smat_add(const struct smat *a, const struct smat *b, struct smat *r)
{
	smat_add(a, b, r);
	debug_op("add", "a", a, "b", b, "r", r);
}

static inline void
__smat_vect(const struct smat *a, const struct smat *v, struct smat *r)
{
	smat_vect(a, v, r);
	debug_op("vect", "a", a, "v", v, "r", r);
}

int main(int argc, char *argv[])
{
	struct smat *a, *b, *r1, *r2, *r3, *v1;
	int i;

	parse_args(argc, argv);
	create_matrices();
	create_vectors();
	create_alphas();
	r1 = smat_calloc(matrix_size, matrix_size);
	r2 = smat_calloc(matrix_size, matrix_size);
	r3 = smat_calloc(matrix_size, matrix_size);
	v1 = smat_calloc(matrix_size, 1);

	a = &matrices[0];
	for (i = 0; i < n_matrices; i++) {
		b = &matrices[i + 1];
		__smat_mult(a, b, r1);
		__smat_scale(a, &alphas[i], r2);
		__smat_add(r1, r2, r3);
		__smat_vect(r3, &vectors[i], v1);
		a = r3;
	}

	/*
	 * SESC does not reliably flush to files, flush here so that
	 * output redirection works on the shell.
	 */
	fflush(NULL);
	return 0;
}
