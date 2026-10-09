# Low-level helpers shared across the package. The flat-index convention here
# is the single source of truth for how (K, D, N) cells map onto the vectors
# passed to Stan; it must match numpy's row-major reshape(-1) of a (K, D, N)
# array, because exposed-cell ordering defines the meaning of treatment_kt[e]
# and the cut-mode provenance columns.

#' Row-major flat index for a (K, D, N) cell
#'
#' Maps 1-based `(k, d, n)` subscripts to the 1-based position the cell
#' occupies in numpy's row-major `reshape(-1)` of a `(K, D, N)` array.
#'
#' @param k,d,n 1-based group / unit / time subscripts (vectorized).
#' @param D,N Number of units and time periods.
#' @return Integer vector of flat indices in `1:(K*D*N)`.
#' @keywords internal
flat_idx <- function(k, d, n, D, N) {
  ((k - 1L) * D + (d - 1L)) * N + n
}

#' Flatten a (K, D, N) R array in row-major (numpy) order
#'
#' R arrays are column-major; this transposes so that element
#' `x[k, d, n]` lands at `flat_idx(k, d, n, D, N)`.
#'
#' @param x A 3-d array with dims (K, D, N).
#' @return A vector of length `K*D*N` in row-major order.
#' @keywords internal
flatten_kdn <- function(x) {
  stopifnot(length(dim(x)) == 3L)
  as.vector(aperm(x, c(3L, 2L, 1L)))
}

#' Inverse of [flatten_kdn()]
#' @param v Vector of length K*D*N in row-major order.
#' @param K,D,N Array dims.
#' @keywords internal
unflatten_kdn <- function(v, K, D, N) {
  aperm(array(v, dim = c(N, D, K)), c(3L, 2L, 1L))
}

#' Recover (k, d, n) subscripts from a row-major flat index
#' @param cell 1-based flat index vector.
#' @param D,N Array dims.
#' @return A list with integer vectors `k`, `d`, `n`.
#' @keywords internal
kdn_from_flat <- function(cell, D, N) {
  c0 <- cell - 1L
  n <- c0 %% N
  rest <- c0 %/% N
  d <- rest %% D
  k <- rest %/% D
  list(k = k + 1L, d = d + 1L, n = n + 1L)
}

#' Numerically stable log-sum-exp
#' @param x Numeric vector.
#' @keywords internal
logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(x - m)))
}

#' Split m selections across chains as evenly as possible
#'
#' Port of Python `cut._chain_quotas`: each chain gets `m %/% n_chains`, and
#' the first `m %% n_chains` chains take one extra.
#'
#' @param m Total number of selections.
#' @param n_chains Number of chains.
#' @return Integer vector of length `n_chains` summing to `m`.
#' @keywords internal
chain_quotas <- function(m, n_chains) {
  base <- m %/% n_chains
  extra <- m %% n_chains
  base + as.integer(seq_len(n_chains) <= extra)
}

#' Evenly strided 1-based indices for output thinning
#'
#' Port of Python `cut.subsample_component_draws`'s index computation:
#' `np.round(np.linspace(0, n - 1, quota)).astype(int)`, returned 1-based.
#' Deterministic (no RNG). numpy rounds half to even, as does R's `round()`.
#'
#' @param n Number of retained draws available.
#' @param quota Number of draws to keep; must satisfy `quota <= n`.
#' @return Sorted integer vector of length `quota` in `1:n`.
#' @keywords internal
strided_indices <- function(n, quota) {
  if (quota > n) {
    cli::cli_abort(
      "Requested {quota} output draw{?s} but only {n} retained draw{?s} available."
    )
  }
  if (quota == 1L) {
    return(1L)
  }
  as.integer(round(seq(0, n - 1, length.out = quota))) + 1L
}

#' Slugify a unit name for filenames
#'
#' Port of Python `tables._slug`: spaces to underscores, lowercased.
#' @param x Character vector.
#' @keywords internal
unit_slug <- function(x) {
  tolower(gsub(" ", "_", x, fixed = TRUE))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Time units the denominator may be weighted by
#'
#' The exposure a rate is quoted per is `denominator * time`, so the time
#' unit is part of the quantity rather than decoration -- changing it changes
#' the number. `"none"` drops the time weighting entirely, which is what a
#' denominator that is itself a per-period flow wants (births during the
#' quarter, say, rather than a population standing through it).
#' @keywords internal
DENOMINATOR_TIME_UNITS <- c("year", "month", "week", "day", "none")

#' Per-row denominator weight, in the requested time unit
#'
#' @param df Rows with `start_date`/`end_date` (else one period each).
#' @param time_unit One of [DENOMINATOR_TIME_UNITS].
#' @keywords internal
time_weight_per_row <- function(df, time_unit = "year") {
  checkmate::assert_choice(time_unit, DENOMINATOR_TIME_UNITS)
  if (identical(time_unit, "none")) {
    return(rep(1, nrow(df)))
  }
  per_year <- switch(time_unit,
    year = 1, month = 12, week = 365.25 / 7, day = 365.25
  )
  years_per_row(df) * per_year
}

#' What one unit of the exposure is called
#'
#' Composed from the denominator's own noun and the time unit:
#' `"person"` + `"year"` gives `"person-years"`, `"birth"` + `"none"` gives
#' `"births"`. The package has no way to know what a denominator counts --
#' population, births, conceptions, vehicle-miles -- so the noun is the
#' user's to supply and defaults to the neutral `"denominator"`. Hardcoding
#' "person-years" baked a demographic assumption into a general-purpose
#' package; dropping the time word entirely hid that rates are per
#' denominator *per unit time*, which is just as misleading.
#'
#' @param denominator_label Singular noun for one unit of the denominator.
#' @param time_unit One of [DENOMINATOR_TIME_UNITS].
#' @keywords internal
format_exposure_noun <- function(denominator_label = "denominator",
                                 time_unit = "year") {
  checkmate::assert_choice(time_unit, DENOMINATOR_TIME_UNITS)
  label <- denominator_label %||% "denominator"
  if (identical(time_unit, "none")) {
    # The field is documented as singular, but "births" is the obvious thing
    # to write for a births denominator and pluralizing it again produced
    # "birthss" in a column header. Take an already-plural noun as given.
    return(if (grepl("s$", label)) label else paste0(label, "s"))
  }
  paste0(label, "-", time_unit, "s")
}

#' Column name for the summed exposure, e.g. `"Person-Years"`
#'
#' Tracks the unit rather than staying fixed: the value under it changes with
#' the time unit, so a stable-but-generic header would mean different things
#' in different runs with nothing in the file to say which.
#' @inheritParams format_exposure_noun
#' @keywords internal
exposure_column_name <- function(denominator_label = "denominator",
                                 time_unit = "year") {
  parts <- strsplit(
    format_exposure_noun(denominator_label, time_unit), "-",
    fixed = TRUE
  )[[1]]
  paste(
    toupper(substring(parts, 1, 1)), substring(parts, 2),
    sep = "", collapse = "-"
  )
}

#' Axis / header text for a rate scaled by `rate_normalizer`
#'
#' "Rate per 1,000 denominator-years", or "Rate per 100,000 birth-years" once
#' the caller says what a unit of the denominator is.
#'
#' @param rate_normalizer Rates are per this many units of exposure.
#' @inheritParams format_exposure_noun
#' @param prefix Leading words, e.g. `"Rate per"` or `"per"`.
#' @param rate_label Complete replacement for the composed text, used
#'   verbatim. Composing from parts keeps the words and the arithmetic in
#'   step, but it cannot know what the numerator counts -- "Rate per 10,000
#'   births" where "deaths per 10,000 births" is what the reader wants -- so
#'   the caller can say it outright instead. Nothing validates it against
#'   `rate_normalizer`, which does change the numbers, so a label naming a
#'   different scale will simply be wrong.
#' @keywords internal
format_rate_label <- function(rate_normalizer,
                              denominator_label = "denominator",
                              time_unit = "year",
                              prefix = "Rate per",
                              rate_label = NULL) {
  if (!is.null(rate_label)) {
    return(rate_label)
  }
  scale <- format(
    rate_normalizer,
    big.mark = ",", scientific = FALSE, trim = TRUE
  )
  paste(
    prefix, scale, format_exposure_noun(denominator_label, time_unit)
  )
}
