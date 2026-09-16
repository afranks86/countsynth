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

#' Axis / header text for a rate scaled by `rate_normalizer`
#'
#' "Rate per 1,000", or "Rate per 1,000 births" when the caller says what a
#' unit of the denominator is. The package has no way to know what a
#' denominator counts -- population, births, conceptions, vehicle-miles -- so
#' the noun is the user's to supply (`output.denominator_label`) and there is
#' none by default. Hardcoding "person-years" baked a demographic assumption
#' into a general-purpose package.
#'
#' @param rate_normalizer Rates are per this many units of exposure.
#' @param denominator_label What one unit of the denominator is, or `NULL`.
#' @param prefix Leading words, e.g. `"Rate per"` or `"per"`.
#' @keywords internal
format_rate_label <- function(rate_normalizer, denominator_label = NULL,
                              prefix = "Rate per") {
  scale <- format(
    rate_normalizer,
    big.mark = ",", scientific = FALSE, trim = TRUE
  )
  paste(c(prefix, scale, denominator_label), collapse = " ")
}
