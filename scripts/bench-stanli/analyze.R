# Summarize out/*.rds: speed, sampling efficiency, and whether the two
# engines' posteriors agree within Monte Carlo error.

`%||%` <- function(a, b) if (is.null(a)) b else a
OUT_DIR <- "/Users/afranks/Dropbox/course/countsynth/scripts/bench-stanli/out"
files <- list.files(OUT_DIR, "\\.rds$", full.names = TRUE)
files <- files[!grepl("^(summary_table|e2e_)", basename(files))]
runs <- lapply(files, readRDS)
# A pair is the two engines' files for one case and seed.
key <- sub("_(cmdstan|stanli)", "", basename(files))

rows <- list()
for (k in unique(key)) {
  pair <- runs[key == k]
  engines <- vapply(pair, `[[`, character(1), "engine")
  if (!setequal(engines, c("cmdstan", "stanli"))) next
  cs <- pair[[which(engines == "cmdstan")]]
  sl <- pair[[which(engines == "stanli")]]
  m <- merge(cs$summary, sl$summary, by = "variable",
             suffixes = c(".cs", ".sl"))
  z <- (m$mean.sl - m$mean.cs) / sqrt(m$mcse_mean.cs^2 + m$mcse_mean.sl^2)
  total <- function(r) r$setup_s + r$sample_s + r$draws_s
  rows[[k]] <- data.frame(
    case = cs$case, seed = if (grepl("_seed", k)) cs$seed else NA,
    cmdstan_s = round(total(cs), 1), stanli_s = round(total(sl), 1),
    speedup = round(total(cs) / total(sl), 2),
    grad_speedup = round(mean(cs$ms_per_grad) / mean(sl$ms_per_grad), 2),
    leap_ratio = round(sum(sl$leap_per_chain) / sum(cs$leap_per_chain), 2),
    med_ess_cs = round(stats::median(m$ess_bulk.cs)),
    med_ess_sl = round(stats::median(m$ess_bulk.sl)),
    min_ess_cs = round(min(m$ess_bulk.cs)), min_ess_sl = round(min(m$ess_bulk.sl)),
    ess_per_s_ratio = round(
      (stats::median(m$ess_bulk.sl) / total(sl)) /
        (stats::median(m$ess_bulk.cs) / total(cs)), 2
    ),
    div_cs = cs$divergent, div_sl = sl$divergent,
    frac_abs_z_gt3 = round(mean(abs(z) > 3, na.rm = TRUE), 4),
    max_abs_z = round(max(abs(z), na.rm = TRUE), 1)
  )
}
tab <- do.call(rbind, rows)
rownames(tab) <- NULL
print(tab[order(tab$case, tab$seed), ], row.names = FALSE)
saveRDS(tab, file.path(OUT_DIR, "summary_table.rds"))

# The null for the z test: the same engine against itself at another seed.
# MCSE-based z-scores over ~1000 correlated cells are not N(0, 1) at 300
# draws, so stanli vs CmdStan is judged against CmdStan vs CmdStan (and
# stanli vs stanli), not against 0.27%.
agree <- function(a, b) {
  m <- merge(a$summary, b$summary, by = "variable")
  z <- (m$mean.x - m$mean.y) / sqrt(m$mcse_mean.x^2 + m$mcse_mean.y^2)
  c(frac_gt3 = mean(abs(z) > 3, na.rm = TRUE), max = max(abs(z), na.rm = TRUE))
}
null_rows <- list()
for (cs_name in unique(vapply(runs, `[[`, character(1), "case"))) {
  of_case <- runs[vapply(runs, function(r) r$case == cs_name, logical(1))]
  if (length(of_case) < 4) next
  pairs <- utils::combn(length(of_case), 2)
  for (j in seq_len(ncol(pairs))) {
    a <- of_case[[pairs[1, j]]]
    b <- of_case[[pairs[2, j]]]
    kind <- paste(sort(c(a$engine, b$engine)), collapse = "-")
    null_rows[[length(null_rows) + 1]] <- data.frame(
      case = cs_name, kind = kind, t(agree(a, b))
    )
  }
}
if (length(null_rows)) {
  nt <- do.call(rbind, null_rows)
  cat("\nPosterior-mean agreement by pair kind (mean over seed pairs):\n")
  print(stats::aggregate(cbind(frac_gt3, max) ~ case + kind, nt, mean),
        digits = 3, row.names = FALSE)
}
