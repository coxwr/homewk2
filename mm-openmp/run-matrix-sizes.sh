#!/bin/sh

setpkgs -a intel_cluster_studio_compiler
CFLAGS="-fast -restrict -std=c99 -openmp"
LDFLAGS="-lm"
# optimal settings, discovered from ./run-chunk-sizes.sh
icc -o matrix_multiply{,.c} $CFLAGS $LDFLAGS -DMAT_SCHEDULE=STATIC_SCHED

outfile=times-matrix-sizes
echo 'size,threads,time' > "$outfile"

export TIMEFORMAT='%E'

for matrix in 2000 4000 8000; do
  echo "$matrix"
  for thread_count in 1 4 8 16 32; do
    echo ">$thread_count"
    export OMP_NUM_THREADS="$thread_count"
    echo -n "$matrix","$thread_count", >> "$outfile"
    (time ./matrix_multiply ./"$matrix".mat "$matrix") 2>>"$outfile" 1>/dev/null
  done
done
