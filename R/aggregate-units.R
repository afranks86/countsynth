# Synthetic aggregate reporting units. Port of aggregate_units.py: appended
# to the draws frame for reporting/PPC only, never fed back into training.
# Per (draw, chain, iteration, time, group): outcome/ypred/denominator are
# summed, treatment is max'd, period boundaries carried, and mu/mu_treated
# pooled by log-sum-exp -- the correct pooling for log-count columns, since
# every estimand downstream sums exp(mu).

source_units_for_spec <- function(spec, source_df) {
  all_units <- unique(source_df$unit)
  missing <- character()
  if (spec$include_all_units) {
    sources <- all_units
  } else if (spec$include_treated_units) {
    sources <- unique(source_df$unit[
      !is.na(source_df$treatment) & source_df$treatment == 1
    ])
  } else {
    requested <- spec$include_units %||% character()
    missing <- setdiff(requested, all_units)
    sources <- intersect(requested, all_units)
  }
  sources <- setdiff(sources, spec$exclude_units %||% character())
  if (length(missing) > 0) {
    if (spec$strict) {
      cli::cli_abort(
        "Aggregate unit {.val {spec$unit}} references missing unit{?s} {.val {missing}}."
      )
    }
    cli::cli_warn(
      "Aggregate unit {.val {spec$unit}} skips missing unit{?s} {.val {missing}}."
    )
  }
  if (length(sources) == 0) {
    cli::cli_warn(
      "Aggregate unit {.val {spec$unit}} resolved to an empty source set; skipping."
    )
  }
  sources
}

logsumexp_narm <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  logsumexp(x)
}

aggregate_one <- function(source_df, unit_name, sources) {
  sub <- source_df[source_df$unit %in% sources, ]
  if (nrow(sub) == 0) {
    return(NULL)
  }
  group_cols <- intersect(
    c(".draw", ".chain", ".iteration", "time", "group"), names(sub)
  )
  has <- function(col) col %in% names(sub)
  # A single source unit's suppressed/missing count would otherwise NA out
  # the whole summed cell (and everything downstream that sums it), even
  # when every other source unit has real data -- e.g. one small state
  # missing a subgroup blanks the entire "all treated states" total for
  # every period it's missing. Fill it with the model's factual expectation
  # instead: mu_treated already equals mu on control cells and mu + te on
  # exposed ones, so exp(mu_treated) is E[Y(1)] where the cell is exposed
  # and E[Y(0)] otherwise -- whichever actually happened there -- and flag
  # the cell as imputed so reporting can mark it. The fill must be constant
  # across draws (the posterior mean, not the per-draw mu_treated): outcome
  # is supposed to be fixed data, and downstream code (compute_quantiles)
  # groups by it to carry it through summaries, assuming it never varies
  # within a (unit, time, group) cell -- a per-draw fill would silently
  # split every draw into its own group.
  outcome_imputed <- FALSE
  if (has("outcome")) {
    outcome_imputed <- is.na(sub$outcome)
    if (has("mu_treated") && any(outcome_imputed)) {
      cell_key <- intersect(c("unit", "time", "group"), names(sub))
      fill <- sub |>
        dplyr::group_by(dplyr::across(dplyr::all_of(cell_key))) |>
        dplyr::summarise(.fill = mean(exp(.data$mu_treated)), .groups = "drop")
      sub <- dplyr::left_join(sub, fill, by = cell_key)
      sub$outcome[outcome_imputed] <- sub$.fill[outcome_imputed]
      sub$.fill <- NULL
    }
  }
  sub$outcome_imputed <- outcome_imputed
  agg <- sub |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::any_of(c("outcome", "ypred", "denominator")), sum
      ),
      dplyr::across(dplyr::any_of("treatment"), max),
      outcome_imputed = any(.data$outcome_imputed),
      dplyr::across(dplyr::any_of(c("start_date", "end_date")), dplyr::first),
      dplyr::across(dplyr::any_of(c("mu", "mu_treated")), logsumexp_narm),
      .groups = "drop"
    )
  agg$unit <- unit_name
  ordered_cols <- intersect(names(source_df), names(agg))
  extra_cols <- setdiff(names(agg), ordered_cols)
  agg[c(ordered_cols, extra_cols)]
}

#' Append synthetic aggregate reporting units to a draws frame
#'
#' Each spec aggregates from the *original* frame, so chaining multiple specs
#' never double-counts. An aggregate colliding with an existing unit is an
#' error unless the spec sets `overwrite`.
#'
#' The names of the units actually created are recorded in an
#' `aggregate_units` attribute on the result, so downstream reporting can tell
#' a pooled unit from a real one -- they are otherwise ordinary rows. Read it
#' with [aggregate_unit_names()].
#'
#' @param draws A `bpnmf_draws` frame.
#' @param specs A list of [bpnmf_aggregate_unit()] specs.
#' @return The draws frame with aggregate-unit rows appended.
#' @export
add_aggregate_units <- function(draws, specs) {
  if (is.null(specs) || length(specs) == 0) {
    return(draws)
  }
  checkmate::assert_list(specs, types = "bpnmf_aggregate_unit")
  attr_names <- c("groups", "units", "times")
  attrs <- stats::setNames(lapply(attr_names, function(a) attr(draws, a)), attr_names)
  attrs <- attrs[!vapply(attrs, is.null, logical(1))]
  source_df <- draws
  result_df <- draws
  existing_units <- unique(draws$unit)
  aggregate_frames <- list()

  for (spec in specs) {
    unit_name <- spec$unit
    if (unit_name %in% existing_units && !spec$overwrite) {
      cli::cli_abort(
        "Aggregate unit {.val {unit_name}} collides with an existing unit;
         set {.field overwrite} = TRUE to allow replacement."
      )
    }
    if (spec$overwrite) {
      result_df <- result_df[result_df$unit != unit_name, ]
      aggregate_frames[[unit_name]] <- NULL
      existing_units <- setdiff(existing_units, unit_name)
    }
    sources <- source_units_for_spec(spec, source_df)
    if (length(sources) == 0) {
      next
    }
    aggregate <- aggregate_one(source_df, unit_name, sources)
    if (is.null(aggregate)) {
      next
    }
    aggregate_frames[[unit_name]] <- aggregate
    existing_units <- c(existing_units, unit_name)
  }

  if (length(aggregate_frames) > 0 && !"outcome_imputed" %in% names(result_df)) {
    result_df$outcome_imputed <- FALSE
  }
  out <- dplyr::bind_rows(c(list(result_df), unname(aggregate_frames)))
  for (nm in names(attrs)) {
    attr(out, nm) <- attrs[[nm]]
  }
  attr(out, "aggregate_units") <- unique(c(
    aggregate_unit_names(draws), names(aggregate_frames)
  ))
  if (!inherits(out, "bpnmf_draws")) {
    class(out) <- c("bpnmf_draws", class(out))
  }
  out
}

#' Names of the synthetic aggregate units in a draws frame
#'
#' Reads the `aggregate_units` attribute stamped by [add_aggregate_units()],
#' which survives row subsetting and `dplyr` verbs. Returns `character(0)` for
#' a frame that never had any.
#'
#' @param draws A `bpnmf_draws` frame.
#' @return Character vector of unit names.
#' @export
aggregate_unit_names <- function(draws) {
  attr(draws, "aggregate_units") %||% character()
}
