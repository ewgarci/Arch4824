//
// Program that tests the performance of the
// Volker Strassen algorithm for matrix multiplication.
//
// W. Cochran  wcochran@vancouver.wsu.edu
//
// Build with optimization:
//  gcc -O3 -std=c99 mm.c -o mm
//
// Two routines compared:
// mmult() : uses convential O(N^3) algorithm.
// mmult_fast() : uses Strassen's O(N^log2(7)) algorithm.
//
// Change preprocessor constant N (defined just above main())
// to a different power of two to try different sizes.
//
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

//define SIZE 37
#define BREAK 1

//
// Classic O(N^3) square matrix multiplication.
// Z = X*Y
// All matrices are NxN and stored in row major order
// each with a specified pitch. 
// The pitch is the distance (in double's) between 
// elements at (row,col) and (row+1,col).
//
void mmult(int N, 
           int Xpitch, const double X[], 
           int Ypitch, const double Y[],
           int Zpitch, double Z[]) {
	int i,j,k;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++) {
      double sum = 0.0;
      for (k = 0; k < N; k++)
        sum += X[i*Xpitch + k]*Y[k*Ypitch + j];
      Z[i*Zpitch + j] = sum;
    }
}

void smat_printf(const char *title, double *M, int N){
	int is_scalar = 0;
	int i, j;

	if (title)
		printf("# name: %s\n", title);
	if (is_scalar) {
		printf("# type: scalar\n");
	} else {
		printf("# type: matrix\n");
		printf("# rows: %u\n", N);
		printf("# columns: %u\n", N);
	}
	
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++)
			printf(" %.17e", M[i*N + j]);
		printf("\n");
	}
}

//
// S = X + Y
//
void madd(int N, 
          int Xpitch, const double X[], 
          int Ypitch, const double Y[],
          int Spitch, double S[]) {
	int i,j;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      S[i*Spitch + j] = X[i*Xpitch + j] + Y[i*Ypitch + j];
}

//
// S = X - Y
//
void msub(int N, 
          int Xpitch, const double X[], 
          int Ypitch, const double Y[],
          int Spitch, double S[]) {
	int i,j;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      S[i*Spitch + j] = X[i*Xpitch + j] - Y[i*Ypitch + j];
}

//
// Volker Strassen algorithm for matrix multiplication.
// Theoretical Runtime is O(N^2.807).
// Assume NxN matrices where N is a power of two.
// Algorithm:
//   Matrices X and Y are split into four smaller
//   (N/2)x(N/2) matrices as follows:
//          _    _          _   _
//     X = | A  B |    Y = | E F |
//         | C  D |        | G H |
//          -    -          -   -
//   Then we build the following 7 matrices (requiring
//   seven (N/2)x(N/2) matrix multiplications -- this is
//   where the 2.807 = log2(7) improvement comes from):
//     P0 = A*(F - H);
//     P1 = (A + B)*H
//     P2 = (C + D)*E
//     P3 = D*(G - E);
//     P4 = (A + D)*(E + H)
//     P5 = (B - D)*(G + H)
//     P6 = (A - C)*(E + F)
//   The final result is
//        _                                            _
//   Z = | (P3 + P4) + (P5 - P1)   P0 + P1              |
//       | P2 + P3                 (P0 + P4) - (P2 + P6)|
//        -                                            -
//
void mmult_fast(int N, 
                int Xpitch, const double X[], 
                int Ypitch, const double Y[],
                int Zpitch, double Z[]) {
				
	int i;
  //
  // Recursive base case.
  // If matrices are 16x16 or smaller we just use
  // the conventional algorithm.
  // At what size we should switch will vary based
  // on hardware platform.
  //
  if (N <= BREAK) {
    mmult(N, Xpitch, X, Ypitch, Y, Zpitch, Z);
    return;
  }

  const int n = N/2;      // size of sub-matrices

  const double *A = X;    // A-D matrices embedded in X
  const double *B = X + n;
  const double *C = X + n*Xpitch;
  const double *D = C + n;

  const double *E = Y;    // E-H matrices embeded in Y
  const double *F = Y + n;
  const double *G = Y + n*Ypitch;
  const double *H = G + n;

  double *P[7];   // allocate temp matrices off heap
  const int sz = n*n*sizeof(double);
  for (i = 0; i < 7; i++)
    P[i] = (double *) malloc(sz);
  double *T = (double *) malloc(sz);
  double *U = (double *) malloc(sz);

  // P0 = A*(F - H);
  msub(n, Ypitch, F, Ypitch, H, n, T);
  mmult_fast(n, Xpitch, A, n, T, n, P[0]);
  
  // P1 = (A + B)*H
  madd(n, Xpitch, A, Xpitch, B, n, T);
  mmult_fast(n, n, T, Ypitch, H, n, P[1]);

  // P2 = (C + D)*E
  madd(n, Xpitch, C, Xpitch, D, n, T);
  mmult_fast(n, n, T, Ypitch, E, n, P[2]);

  // P3 = D*(G - E);
  msub(n, Ypitch, G, Ypitch, E, n, T);
  mmult_fast(n, Xpitch, D, n, T, n, P[3]);

  // P4 = (A + D)*(E + H)
  madd(n, Xpitch, A, Xpitch, D, n, T);
  madd(n, Ypitch, E, Ypitch, H, n, U);
  mmult_fast(n, n, T, n, U, n, P[4]);

  // P5 = (B - D)*(G + H)
  msub(n, Xpitch, B, Xpitch, D, n, T);
  madd(n, Ypitch, G, Ypitch, H, n, U);
  mmult_fast(n, n, T, n, U, n, P[5]);

  // P6 = (A - C)*(E + F)
  msub(n, Xpitch, A, Xpitch, C, n, T);
  madd(n, Ypitch, E, Ypitch, F, n, U);
  mmult_fast(n, n, T, n, U, n, P[6]);

  // Z upper left = (P3 + P4) + (P5 - P1)
  madd(n, n, P[4], n, P[3], n, T);
  msub(n, n, P[5], n, P[1], n, U);
  madd(n, n, T, n, U, Zpitch, Z);

  // Z lower left = P2 + P3
  madd(n, n, P[2], n, P[3], Zpitch, Z + n*Zpitch);

  // Z upper right = P0 + P1
  madd(n, n, P[0], n, P[1], Zpitch, Z + n);
  
  // Z lower right = (P0 + P4) - (P2 + P6)
  madd(n, n, P[0], n, P[4], n, T);
  madd(n, n, P[2], n, P[6], n, U);
  msub(n, n, T, n, U, Zpitch, Z + n*(Zpitch + 1));

  free(U);  // deallocate temp matrices
  free(T);
  for (i = 6; i >= 0; i--)
    free(P[i]);
}

void mprint(int N, int pitch, const double M[]) {

	int i,j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
      printf("%+0.4f ", M[i*pitch + j]);
    printf("\n");
  }
}

void mrand(int N, int pitch, double M[]) {
  //const double r = 10.0;
  int i,j;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      M[i*pitch + j] = random()/(RAND_MAX * 2.0);
}


/*
int timeval_subtract(struct timeval *result, 
                     struct timeval *t2, struct timeval *t1) {
  long int diff = 
    (t2->tv_usec + 1000000 * t2->tv_sec) - 
    (t1->tv_usec + 1000000 * t1->tv_sec);
  result->tv_sec = diff / 1000000;
  result->tv_usec = diff % 1000000;
  return (diff<0);
}

void timeval_print(struct timeval *tv) {
  char buffer[30];
  time_t curtime;

  printf("%ld.%06ld", (long int) tv->tv_sec, (long int) tv->tv_usec);
  curtime = tv->tv_sec;
  strftime(buffer, 30, "%m-%d-%Y  %T", localtime(&curtime));
  printf(" = %s.%06ld\n", buffer, (long int) tv->tv_usec);
}
*/




int main(int argc, char *argv[]) {

  double *X, *Y, *Z, *Zfast;
  int SIZE = atoi(argv[1]);
  X = (double*) malloc(SIZE*SIZE*sizeof(double));
  Y = (double*) malloc(SIZE*SIZE*sizeof(double));
  Z = (double*) malloc(SIZE*SIZE*sizeof(double));
  Zfast = (double*) malloc(SIZE*SIZE*sizeof(double));

  mrand(SIZE, SIZE, X);
  mrand(SIZE, SIZE, Y);
  mrand(SIZE, SIZE, Zfast);


  mmult_fast(SIZE, SIZE, X, SIZE, Y, SIZE, Zfast);

	printf("# op: mult\n");
	smat_printf("a", X, SIZE);
	smat_printf("b", Y, SIZE);
	smat_printf("r", Zfast, SIZE);
	printf("# ---\n");
	
  return 0;
}
