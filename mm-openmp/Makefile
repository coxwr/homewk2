CC := setpkgs -a intel_cluster_studio_compiler && icc

CFLAGS += -fast -restrict -std=c99 -openmp
LDFLAGS += -lm

matrix_multiply : matrix_multiply.c
	$(CC) -o $@ $< $(CFLAGS) $(LDFLAGS)

generate_rand_matrix : generate_rand_matrix.c
	$(CC) -o $@ $< $(CFLAGS)

.PHONY: clean

clean :
	rm matrix_multiply
