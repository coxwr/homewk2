//*****************************************************************************
//
// This program calculates the product of a square matrix with itself:
//
// B = A * A
//
// Please keep all code in this single file.
//
//
//*****************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

const char *HELP_TEXT =
    "This program computes the product of an n x n matrix with itself\n"
    "Usage: ./matrix_multiply FILE N";
const size_t BASE = 10;

static inline int index_matrix(const int *restrict mat, size_t x, size_t y,
                               size_t n) {
  return mat[x * n + y];
}

static inline int *index_matrix_set(int *restrict mat, size_t x, size_t y,
                                    size_t n) {
  return mat + (x * n + y);
}

int parse_from_file(const char *path, int *restrict mat, size_t n) {
  /* assume n > 0 */
  if (!path || !mat) {
    return 1;
  }
  FILE *input = fopen(path, "r");
  if (!input) {
    return 1;
  }
  int retc = 0;
  for (size_t i = 0; i < n * n; ++i) {
    if (fscanf(input, "%d", mat + i) != 1) {
      retc = 1;
      goto cleanup;
    }
  }
cleanup:
  if (fclose(input)) {
    retc = 1;
  }
  return retc;
}

static inline int dot_product_row_col(int *restrict a, size_t n, size_t row,
                                      size_t col) {
  int sum = 0;
  for (size_t i = 0; i < n; ++i) {
    sum += index_matrix(a, row, i, n) * index_matrix(a, i, col, n);
  }
  return sum;
}

#define STATIC_SCHED 1
#define DYNAMIC_SCHED 2
#define GUIDED_SCHED 3
#ifndef MAT_CHUNK_SIZE
#define MAT_CHUNK_ARG
#else
#define MAT_CHUNK_ARG , MAT_CHUNK_SIZE
#endif
void square_matrix_into(int *restrict a, int *restrict b, size_t n) {
#if MAT_SCHEDULE == STATIC_SCHED
#pragma omp parallel for schedule(static MAT_CHUNK_ARG)
#elif MAT_SCHEDULE == DYNAMIC_SCHED
#pragma omp parallel for schedule(dynamic MAT_CHUNK_ARG)
#elif MAT_SCHEDULE == GUIDED_SCHED
#pragma omp parallel for schedule(guided MAT_CHUNK_ARG)
#else
#error "no schedule type given!"
#endif
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      *index_matrix_set(b, i, j, n) = dot_product_row_col(a, n, i, j);
    }
  }
}

int write_matrix_stdout(int *restrict mat, size_t n) {
  for (size_t i = 0; i < n * n; ++i) {
    if (printf("%d ", mat[i]) == 0) {
      return 1;
    }
    if ((i + 1) % n == 0 && printf("\n") == 0) {
      return 1;
    }
  }
  return 0;
}

int main(int argc, char **argv) {
  int retc = 0;
  // check command line arguments
  if (argc != 3) {
    fprintf(stderr, "%s\n", HELP_TEXT);
    retc = 1;
    goto exit;
  }

  // TODO: parse input arguments
  FILE *input = fopen(argv[1], "r");
  if (!input) {
    fprintf(stderr, "%s is not a valid filename\n", argv[1]);
    retc = 1;
    goto exit;
  }
  char *tailptr;
  long n = strtol(argv[2], &tailptr, BASE);
  /* strtol failed, or found negative */
  if (n <= 0 || (0 == n && tailptr == argv[2])) {
    fprintf(stderr, "invalid N given: %s\n", argv[2]);
    retc = 1;
    goto exit;
  }

  // TODO: dynamically allocate space for matrix_A (input matrix) in 1d array
  int *matrix_A = malloc(n * n * sizeof(int)); // declare input matrix
  // TODO: dynamically allocate space for matrix_B (output matrix) in 1d array
  int *matrix_B = malloc(n * n * sizeof(int)); // declare output matrix

  // TODO: call function to read data from file and copy into matrix_A
  if (parse_from_file(argv[1], matrix_A, n)) {
    fprintf(stderr, "file parsing failed from %s\n", argv[1]);
    retc = 1;
    goto cleanup;
  }

  // TODO: call function to perform matrix multiplication ( matrix_B = matrix_A
  // * matrix_A )
  square_matrix_into(matrix_A, matrix_B, n);

  // TODO: call function to write results (matrix_B) to stdout
  if (write_matrix_stdout(matrix_B, n)) {
    fprintf(stderr, "writing matrix to stdout failed\n");
    retc = 1;
    goto cleanup;
  }

// TODO: free space allocated for matrix_A and matrix_B
cleanup:
  if (fclose(input)) {
    perror("input file could not be closed");
  }
  free(matrix_A);
  free(matrix_B);
exit:
  return retc;
}
