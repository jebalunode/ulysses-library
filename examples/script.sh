awk '
/cycles/ {gsub(",","",$1); cycles=$1}
/instructions/ {gsub(",","",$1); instr=$1}
/fp_ret_sse_avx_ops.all/ {gsub(",","",$1); fp=$1}
END {
  if (cycles > 0 && instr > 0) {
    ipc = instr / cycles
    vec = (fp / instr) * 100
    flops = fp * 8
    print "IPC:", ipc
    print "Vectorization %:", vec
    print "Approx FLOPS (assuming AVX2 8 floats):", flops
  } else {
    print "Error: missing perf counters"
  }
}'
