#!/bin/sh

setpkgs -a intel_cluster_studio_compiler
export TIMEFORMAT='%E'

outfile=sim-trials
echo "size,threads,time" > "$outfile"

for size in 125 1000 4000; do
  for threads in 1 4 8 16 32; do
    make clean && make
    echo -n "$size","$threads", >> "$outfile"
    export OMP_NUM_THREADS="$threads"
    (time bin/run_md -N "$size") 2>>"$outfile" 1>/dev/null
  done
done
