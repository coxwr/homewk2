#!/bin/sh

setpkgs -a intel_cluster_studio_compiler
CFLAGS="-fast -restrict -std=c99 -openmp"
LDFLAGS="-lm"

export TIMEFORMAT='%E'

outfile=times-chunk-sizes
echo 'sched,chunk,time' > "$outfile"

for sched_type in STATIC_SCHED DYNAMIC_SCHED GUIDED_SCHED; do
  echo "$sched_type"
  echo ">-1"
  icc -o matrix_multiply{,.c} $CFLAGS $LDFLAGS -DMAT_SCHEDULE="$sched_type"
  echo -n "$sched_type","-1", >>"$outfile"
  (time ./matrix_multiply ./4000.mat 4000) 2>>"$outfile" 1>/dev/null
  for chunk_size in 5 50 500 1000 2000; do
    echo ">$chunk_size"
    icc -o matrix_multiply{,.c} $CFLAGS $LDFLAGS -DMAT_SCHEDULE="$sched_type" \
        -DMAT_CHUNK_SIZE="$chunk_size"
    echo -n "$sched_type","$chunk_size", >>"$outfile"
    (time ./matrix_multiply ./4000.mat 4000) 2>>"$outfile" 1>/dev/null
  done
done
