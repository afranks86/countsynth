# Synthetic aggregate reporting units. Port of aggregate_units.py: appended
# to the draws frame for reporting/PPC only, never fed back into training.
# Per (draw, chain, iteration, time, group): outcome/ypred/denominator are
# summed, treatment is max'd, period boundaries carried, and mu/mu_treated
# pooled by log-sum-exp — the correct pooling for log-count columns, since
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
  agg <- sub |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::any_of(c("outcome", "ypred", "denominator")), sum
      ),
      dplyr::across(dplyr::any_of("treatment"), max),
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

  out <- dplyr::bind_rows(c(list(result_df), unname(aggregate_frames)))
  for (nm in names(attrs)) {
    attr(out, nm) <- attrs[[nm]]
  }
  if (!inherits(out, "bpnmf_draws")) {
    class(out) <- c("bpnmf_draws", class(out))
  }
  out
}
