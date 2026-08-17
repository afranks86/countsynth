# Bundled example: the US fertility panel (state x bimonth births by
# subgroup) shipped with the Python package, with a smoke-test configuration
# mirroring configs/fertility_smoke_test.yaml.

#' Path to the bundled fertility example data
#' @export
bpnmf_example_data <- function() {
  system.file("extdata", "fertility_data.csv", package = "bpnmf")
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
#' @param ... Overrides passed to [bpnmf_config()] fields (e.g. `mcmc =`).
#' @return A [bpnmf_config()] object.
#' @export
bpnmf_example_config <- function(output_dir = file.path(tempdir(), "bpnmf_results"),
                                 ...) {
  args <- list(
    input_file = bpnmf_example_data(),
    output_dir = output_dir,
    schema = bpnmf_schema(
      unit_col = "state", time_col = "time", treatment_col = "exposed_births",
      outcomes_from_prefixes = bpnmf_prefixes(
        outcome_prefix = "births_", denominator_prefix = "pop_",
        include = "total"
      )
    ),
    model = bpnmf_model_opts(
      outcome_distribution = "NB",
      types = list(
        total = bpnmf_type(groups = "total", ranks_to_test = 3)
      )
    ),
    mcmc = bpnmf_mcmc_opts(
      iter_warmup = 200, iter_sampling = 200, thin = 1,
      seed = 8675309, progress = FALSE
    ),
    output = bpnmf_output_opts(figures = TRUE),
    date_format = "%Y-%m-%d",
    start_date = "2016-01-01", end_date = "2024-01-01",
    aggregation = bpnmf_aggregation(enabled = TRUE, period = "bimonthly")
  )
  override <- list(...)
  for (nm in names(override)) {
    args[[nm]] <- override[[nm]]
  }
  do.call(bpnmf_config, args)
}
