#!/bin/sh
# Package-level runs at each config's own settings, backends interleaved.
# /usr/bin/time -l reports the R process's peak RSS (CmdStan's chains are
# separate processes and not counted; stanli's are threads and are).
cd "$(dirname "$0")"
run() {
  /usr/bin/time -l Rscript e2e.R "$@" > out/e2e_tmp.log 2>&1
  grep -E "^E2E|Error|rror in" out/e2e_tmp.log
  grep "maximum resident set size" out/e2e_tmp.log | awk -v b="$4" -v c="$1" '{printf "RSS %s %s %.0f MB\n", c, b, $1/1048576}'
}
for e in cmdstan stanli; do run configs/mortality_full.yml race 5 $e; done
for e in cmdstan stanli; do run configs/mortality_con.yml con_noncon 3 $e; done
for e in cmdstan stanli; do run configs/fertility_te_covariates.yml total 6 $e; done
echo DONE
