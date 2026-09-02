# Shared test setup.

has_cmdstan <- tryCatch(
  nzchar(cmdstanr::cmdstan_path()),
  error = function(e) FALSE
)

skip_if_no_cmdstan <- function() {
  testthat::skip_if_not(has_cmdstan, "CmdStan not available")
}

fixture_path <- function(name) {
  testthat::test_path("fixtures", name)
}

# Tiny synthetic panel shared across tests: K=2 groups, D=3 units, N=5 times,
# unit "C" treated from period 4 on.
make_test_data <- function(seed = 42) {
  withr::with_seed(seed, {
    groups <- c("g1", "g2")
    units <- c("A", "B", "C")
    times <- seq(as.Date("2020-01-01"), by = "month", length.out = 5)
    grid <- expand.grid(
      unit = units, time = times, group = groups,
      stringsAsFactors = FALSE
    )
    grid$treatment <- as.integer(grid$unit == "C" & grid$time >= times[4])
    grid$denominator <- runif(nrow(grid), 5000, 30000)
    grid$outcome <- rpois(nrow(grid), 50 * grid$denominator / 1e4)
    tibble::as_tibble(grid)
  })
}

make_test_config <- function(df_path, output_dir,
                             inference_mode = NULL, ...) {
  bpnmf_config(
    input_file = df_path,
    output_dir = output_dir,
    schema = bpnmf_schema(
      unit_col = "unit", time_col = "time", treatment_col = "treatment",
      outcomes = list(
        bpnmf_outcome("outcome_g1", "g1", denominator_col = "denom_g1"),
        bpnmf_outcome("outcome_g2", "g2", denominator_col = "denom_g2")
      )
    ),
    model = bpnmf_model_opts(
      outcome_distribution = "NB",
      types = list(
        both = bpnmf_type(groups = c("g1", "g2"), ranks_to_test = 2)
      ),
      inference_mode = inference_mode
    ),
    mcmc = bpnmf_mcmc_opts(
      auto_parallelism = FALSE, chains = 2, parallel_chains = 2,
      iter_warmup = 150, iter_sampling = 150, thin = 1, seed = 99
    ),
    ...
  )
}

# Write the long test panel as a wide CSV matching make_test_config's schema.
write_test_csv <- function(df, path) {
  wide <- df |>
    tidyr::pivot_wider(
      id_cols = c("unit", "time", "treatment"),
      names_from = "group",
      values_from = c("outcome", "denominator"),
      names_glue = "{ifelse(.value == 'outcome', 'outcome', 'denom')}_{group}"
    )
  utils::write.csv(wide, path, row.names = FALSE)
  path
}

# 2 units (A control, B treated in periods 2-3), 1 group, 3 times, 4 draws.
make_draws_frame <- function() {
  # 2 units (A control, B treated in periods 2-3), 1 group, 3 times, 4 draws
  units <- c("A", "B")
  times <- seq(as.Date("2021-01-01"), by = "month", length.out = 3)
  grid <- expand.grid(
    unit = units, time = times,
    stringsAsFactors = FALSE
  )
  grid$group <- "total"
  grid$treatment <- as.integer(grid$unit == "B" & grid$time >= times[2])
  grid$outcome <- c(100, 110, 100, 120, 100, 130)
  grid$denominator <- 1000
  draws <- dplyr::bind_rows(lapply(1:4, function(d) {
    g <- grid
    g$.draw <- d
    g$.chain <- 1L
    g$.iteration <- d
    g$mu <- log(100 + d)             # untreated log-count
    g$mu_treated <- log(100 + d) + ifelse(g$treatment == 1, log(1.2), 0)
    g$ypred <- 100 + d
    g
  }))
  class(draws) <- c("bpnmf_draws", class(draws))
  draws
}
