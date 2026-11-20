#!/bin/sh
# Engines interleaved per case so drift (thermal, background load) hits both.
cd "$(dirname "$0")"
for c in fertility_full mortality_con mortality_full fertility_te_race; do
  for e in cmdstan stanli; do
    Rscript sample.R "$c" "$e" "${1:-300}" "${2:-300}" 2>&1 | grep -E "^(fertility|mortality)|Error"
  done
done
echo DONE
