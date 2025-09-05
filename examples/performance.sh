#!/bin/bash
perf stat -e cycles,instructions,fp_ret_sse_avx_ops.all "$@" >& log1

cat log1 | \
awk '
/cycles/ {gsub(",","",$1); cycles=$1}
/instructions/ {gsub(",","",$1); instr=$1}
/fp_ret_sse_avx_ops.all/ {gsub(",","",$1); fp=$1}
/elapsed/ {gsub(",","",$1); time=$1}
END {
  if (cycles > 0 && instr > 0) {
    ipc = instr / cycles
    vec = (fp / instr) * 100
    flopsf = fp * 8
    flopsd = fp *4
    print "IPC:", ipc
    print "Vectorization %:", vec
    print "Total runtime %:", time
    print "Approx FLOPS (assuming AVX2 8 floats):", flopsf/time
    print "Approx GigaFLOPS (assuming AVX2 8 floats):", flopsf/1000000000.0/time
    print "Approx FLOPS (assuming AVX2 4 doubles):", flopsd/time
    print "Approx GigaFLOPS (assuming AVX2 4 doubles):", flopsd/1000000000.0/time


  } else {
    print "Error: missing perf counters"
  }
}'

