#!/bin/sh
# Replicate seeds: is stanli's lower min-ESS in the first pass noise, or does
# its NUTS mix differently? Engines interleaved per seed.
cd "$(dirname "$0")"
for s in 101 202 303; do
  for c in fertility_full mortality_con mortality_full; do
    for e in cmdstan stanli; do
      Rscript sample.R "$c" "$e" 300 300 "$s" 2>&1 | grep -E "^(fertility|mortality)|Error"
    done
  done
done
echo DONE
