# Bundled example: the US fertility panel (state x bimonth births by
# subgroup) shipped with the Python package, with a smoke-test configuration
# mirroring configs/fertility_smoke_test.yaml.

#' Path to the bundled fertility example data
#' @export
countsynth_example_data <- function() {
  system.file("extdata", "fertility_data.csv", package = "countsynth")
}

#' Smoke-test configuration for the bundled fertility data
#'
#' Mirrors the Python package's `fertility_smoke_test.yaml`: total births,
#' rank 3, bimonthly aggregation, 200 warmup + 200 sampling iterations.
#' Fast enough to verify the pipeline end-to-end (~minutes); not for
#' inference.
#'
#' @param output_dir Where run artifacts are written (default: a session
#'   temporary directory).
#' @param ... Overrides passed to [countsynth_config()] fields (e.g. `mcmc =`).
#' @return A [countsynth_config()] object.
#' @export
countsynth_example_config <- function(output_dir = file.path(tempdir(), "countsynth_results"),
                                 ...) {
  args <- list(
    input_file = countsynth_example_data(),
    output_dir = output_dir,
    schema = countsynth_schema(
      unit_col = "state", time_col = "time", treatment_col = "exposed_births",
      outcomes_from_prefixes = countsynth_prefixes(
        outcome_prefix = "births_", denominator_prefix = "pop_",
        include = "total"
      )
    ),
    model = countsynth_model_opts(
      outcome_distribution = "NB",
      types = list(
        total = countsynth_type(groups = "total", ranks_to_test = 3)
      )
    ),
    mcmc = countsynth_mcmc_opts(
      iter_warmup = 200, iter_sampling = 200, thin = 1,
      seed = 8675309, progress = FALSE
    ),
    output = countsynth_output_opts(figures = TRUE),
    date_format = "%Y-%m-%d",
    start_date = "2016-01-01", end_date = "2024-01-01",
    time_aggregation = countsynth_time_aggregation(enabled = TRUE, period = "bimonthly")
  )
  override <- list(...)
  for (nm in names(override)) {
    args[[nm]] <- override[[nm]]
  }
  do.call(countsynth_config, args)
}
