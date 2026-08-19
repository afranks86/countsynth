# Treatment-effect regression: formula parsing, design construction, and the
# Stan data block. The design is the contract between the R formula and the
# Stan `te` assembly, so these tests pin the row ordering (canonical exposed
# flat order), the ragged random-effect layout, and the legacy-path default.

# A panel with two treated units adopting at different periods: unit-level
# covariates then vary across the exposed cells, and event_time is staggered.
# K=2 groups, D=4 units, N=6 times; C treated from period 4, D from period 5.
make_te_data <- function() {
  groups <- c("g1", "g2")
  units <- c("A", "B", "C", "D")
  times <- seq(as.Date("2020-01-01"), by = "month", length.out = 6)
  grid <- expand.grid(
    unit = units, time = times, group = groups, stringsAsFactors = FALSE
  )
  grid$treatment <- as.integer(
    (grid$unit == "C" & grid$time >= times[4]) |
      (grid$unit == "D" & grid$time >= times[5])
  )
  grid$denominator <- 1e4
  grid$outcome <- 50
  arrays <- build_model_arrays(grid, groups)
  arrays$df <- grid
  arrays$type <- "test"
  structure(arrays, class = c("bpnmf_data", "list"))
}

test_that("covariate frame has one row per exposed cell in flat order", {
  data <- make_te_data()
  frame <- te_covariate_frame(data)
  exp_cell <- which(!flatten_kdn(data$control_idx_array))
  sub <- kdn_from_flat(exp_cell, length(data$units), length(data$times))

  expect_equal(nrow(frame), length(exp_cell))
  expect_equal(as.character(frame$group), data$groups[sub$k])
  expect_equal(as.character(frame$unit), data$units[sub$d])
  expect_equal(frame$time, data$times[sub$n])
  expect_equal(frame$time_idx, sub$n)
  # Staggered adoption: C is treated from period 4 (event_time 0..2), D from
  # period 5 (0..1), so event_time is per-unit, not calendar time.
  expect_equal(sort(unique(frame$event_time)), c(0, 1, 2))
  expect_equal(sort(unique(frame$event_time[frame$unit == "D"])), c(0, 1))
  expect_equal(
    frame$event_time[frame$unit == "C" & frame$time_idx == 4][[1]], 0
  )
  expect_equal(
    frame$event_time[frame$unit == "D" & frame$time_idx == 6][[1]], 1
  )
})

test_that("fixed-effect design matches model.matrix on the same variables", {
  data <- make_te_data()
  design <- build_te_design(
    bpnmf_te_opts(formula = ~ 1 + event_time, standardize = FALSE), data
  )
  frame <- te_covariate_frame(data)
  expect_equal(design$x_names, c("(Intercept)", "event_time"))
  expect_equal(unname(design$X[, 2]), as.numeric(frame$event_time))
  expect_equal(length(design$terms), 0L)
})

test_that("random-effect terms carry levels, indices, and their own design", {
  data <- make_te_data()
  design <- build_te_design(
    bpnmf_te_opts(
      formula = ~ 1 + event_time + (1 + event_time | group),
      standardize = FALSE
    ),
    data
  )
  expect_equal(length(design$terms), 1L)
  term <- design$terms[[1]]
  expect_equal(term$label, "group")
  expect_equal(term$levels, data$groups)
  expect_equal(term$predictors, c("(Intercept)", "event_time"))
  frame <- te_covariate_frame(data)
  expect_equal(term$index, as.integer(frame$group))
  expect_equal(unname(term$Z[, 2]), as.numeric(frame$event_time))
})

test_that("interaction grouping factors are supported", {
  data <- make_te_data()
  design <- build_te_design(
    bpnmf_te_opts(formula = ~ 1 + (1 | group:unit)), data
  )
  term <- design$terms[[1]]
  expect_equal(term$label, "group:unit")
  # Only the (group, unit) combinations that actually appear are levels.
  frame <- te_covariate_frame(data)
  present <- unique(paste(frame$group, frame$unit, sep = ":"))
  expect_setequal(term$levels, present)
})

test_that("standardization centers and scales continuous covariates only", {
  data <- make_te_data()
  covs <- data.frame(unit = data$units, dist = c(0, 10, 20, 30))
  design <- build_te_design(
    bpnmf_te_opts(formula = ~ 1 + dist, covariates = covs), data
  )
  expect_named(design$centers, "dist")
  expect_equal(mean(design$frame$dist), 0, tolerance = 1e-12)
  expect_equal(stats::sd(design$frame$dist), 1, tolerance = 1e-12)
  # raw_frame keeps the original scale so plots can map back exactly.
  expect_equal(
    design$frame$dist * design$scales[["dist"]] + design$centers[["dist"]],
    design$raw_frame$dist
  )
  # Only the exposed units (C, D) appear.
  expect_setequal(design$raw_frame$dist, c(20, 30))
})

test_that("indicator covariates keep their 0/1 coding", {
  data <- make_te_data()
  covs <- data.frame(unit = data$units, flag = c(0, 0, 1, 0))
  design <- build_te_design(
    bpnmf_te_opts(formula = ~ 1 + flag, covariates = covs), data
  )
  expect_false("flag" %in% names(design$centers))
  expect_setequal(design$frame$flag, c(0, 1))
})

test_that("a rank-deficient fixed design is refused", {
  data <- make_te_data()
  # Constant across the exposed cells: only unit C and D are exposed and both
  # carry the same value, so the column duplicates the intercept.
  covs <- data.frame(unit = data$units, dist = c(0, 10, 20, 20))
  expect_error(
    build_te_design(
      bpnmf_te_opts(formula = ~ 1 + dist, covariates = covs), data
    ),
    "rank deficient"
  )
  # Collinear pair.
  covs2 <- data.frame(
    unit = data$units, a = c(0, 1, 2, 3), b = c(0, 2, 4, 6)
  )
  expect_error(
    build_te_design(
      bpnmf_te_opts(formula = ~ 0 + a + b, covariates = covs2), data
    ),
    "rank deficient"
  )
})

test_that("covariates join on whichever key columns they carry", {
  data <- make_te_data()
  by_unit_time <- expand.grid(
    unit = data$units, time = data$times, stringsAsFactors = FALSE
  )
  by_unit_time$w <- seq_len(nrow(by_unit_time))
  frame <- te_covariate_frame(data, by_unit_time)
  expect_false(anyNA(frame$w))
  keyed <- paste(frame$unit, frame$time)
  cov_key <- paste(by_unit_time$unit, by_unit_time$time)
  ref <- by_unit_time$w[match(keyed, cov_key)]
  expect_equal(frame$w, ref)
})

test_that("covariate problems are caught with actionable messages", {
  data <- make_te_data()
  expect_error(
    build_te_design(bpnmf_te_opts(formula = ~nosuchvar), data),
    "neither a built-in variable"
  )
  expect_error(
    build_te_design(bpnmf_te_opts(formula = ~ 1 + (1 | event_time)), data),
    "must be categorical"
  )
  expect_error(
    build_te_design(bpnmf_te_opts(formula = ~time), data),
    "Date column"
  )
  # Missing on some exposed cells: covariate defined for one unit only.
  partial <- data.frame(unit = data$units[[1]], z = 1)
  expect_error(
    build_te_design(
      bpnmf_te_opts(formula = ~z, covariates = partial), data
    ),
    "is missing on"
  )
  expect_error(
    bpnmf_te_opts(covariates = data.frame(unit = "A", z = 1)),
    "without a .*formula"
  )
  expect_error(
    build_te_design(
      bpnmf_te_opts(
        formula = ~z,
        covariates = data.frame(nokey = "A", z = 1)
      ),
      data
    ),
    "at least one key column"
  )
})

test_that("te_stan_fields lays out the ragged random-effect blocks", {
  data <- make_te_data()
  design <- build_te_design(
    bpnmf_te_opts(
      formula = ~ 1 + event_time + (1 + event_time | group) + (1 | unit),
      standardize = FALSE
    ),
    data
  )
  sd <- te_stan_fields(design, design$n_exposed)
  # Only exposed units can carry a treatment-effect random intercept, so the
  # (1 | unit) term has one level per TREATED unit, not per panel unit.
  n_treated <- length(unique(te_covariate_frame(data)$unit))
  expect_equal(n_treated, 2L)
  expect_equal(sd$te_reg, 1L)
  expect_equal(sd$P, 2L)
  expect_equal(sd$J, 2L)
  expect_equal(as.integer(sd$Q), c(2L, 1L))
  expect_equal(as.integer(sd$L), c(length(data$groups), n_treated))
  expect_equal(sd$Qtot, 3L)
  expect_equal(sd$Utot, 2L * length(data$groups) + n_treated)
  # Z is the terms' designs concatenated by column; re_level one row per term.
  expect_equal(ncol(sd$Z), sd$Qtot)
  expect_equal(dim(sd$re_level), c(2L, design$n_exposed))
  expect_equal(length(sd$te_beta_prior_scale), sd$P)
  expect_equal(length(sd$te_re_prior_scale), sd$Qtot)
})

test_that("no design gives the legacy stan-data block", {
  data <- make_te_data()
  n_exposed <- sum(!data$control_idx_array)
  sd <- te_stan_fields(NULL, n_exposed)
  expect_equal(sd$te_reg, 0L)
  expect_equal(sd$P, 0L)
  expect_equal(sd$J, 0L)
  expect_equal(sd$Utot, 0L)

  # And the full stan-data builders default to it.
  joint <- stan_data_joint(data, rank = 2, model_treated = TRUE)
  expect_equal(joint$te_reg, 0L)
  stage2 <- stan_data_stage2(
    data,
    mu_ctrl_flat = rep(0, length(flatten_kdn(data$Y))),
    phi_unit = rep(1, length(data$units))
  )
  expect_equal(stage2$te_reg, 0L)
  expect_equal(nrow(stage2$X), stage2$n_exposed)
})

test_that("a design requires model_treated", {
  data <- make_te_data()
  design <- build_te_design(bpnmf_te_opts(formula = ~ 1 + event_time), data)
  expect_error(
    stan_data_joint(data, rank = 2, model_treated = FALSE, te_design = design),
    "requires .*model_treated"
  )
  expect_error(
    bpnmf_model_opts(
      model_treated = FALSE,
      treatment_effects = bpnmf_te_opts(formula = ~event_time)
    ),
    "requires .*model_treated"
  )
})

test_that("the YAML loader round-trips a treatment_effects section", {
  tmp <- withr::local_tempdir()
  yml <- file.path(tmp, "config.yaml")
  writeLines(c(
    "data:",
    "  input_file: in.csv",
    "  output_dir: out",
    "  schema:",
    "    unit_col: unit",
    "    time_col: time",
    "    treatment_col: treatment",
    "    outcomes:",
    "      - {outcome_col: outcome_g1, label: g1}",
    "model:",
    "  types:",
    "    both: {groups: [g1], ranks_to_test: [2]}",
    "  treatment_effects:",
    "    formula: ~ 1 + event_time + (1 + event_time | group)",
    "    standardize: false",
    "    coef_prior_scale: 2.5"
  ), yml)
  config <- read_bpnmf_config(yml)
  te <- config$model$treatment_effects
  expect_s3_class(te, "bpnmf_te_opts")
  expect_false(te$standardize)
  expect_equal(te$coef_prior_scale, 2.5)
  expect_equal(
    deparse1(te$formula),
    deparse1(~ 1 + event_time + (1 + event_time | group))
  )
})

test_that("unknown treatment_effects keys are rejected", {
  tmp <- withr::local_tempdir()
  yml <- file.path(tmp, "config.yaml")
  writeLines(c(
    "data:",
    "  input_file: in.csv",
    "  output_dir: out",
    "  schema:",
    "    unit_col: unit",
    "    time_col: time",
    "    treatment_col: treatment",
    "    outcomes:",
    "      - {outcome_col: outcome_g1, label: g1}",
    "model:",
    "  treatment_effects:",
    "    formula: ~ event_time",
    "    bogus_key: 1"
  ), yml)
  expect_error(read_bpnmf_config(yml), "bogus_key")
})
