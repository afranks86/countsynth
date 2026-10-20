# Synthetic aggregate reporting groups: the group-axis counterpart to
# add_aggregate_units(). Pooling races/ages/causes into an implied total is
# reporting-only, exactly like a pooled unit, and shares aggregate_over().
#
# This is NOT the same as putting "total" in a model type's `groups` (see
# total_from / total_all in the data schema). That fits the total as another
# modeled group: it re-enters the likelihood, double-counts the same events,
# and consumes a factor dimension. A pooled reporting group never touches the
# fit -- it reads the posterior that is already there.
#
# The pooled effect it reports is sum(exp(mu_treated)) / sum(exp(mu)) over the
# member groups, i.e. the total weighted by each group's own exposure. That is
# the implied overall effect, not the average of the per-group percentages --
# a large group moves it more than a small one, as it should.

source_groups_for_spec <- function(spec, source_df) {
  all_groups <- unique(source_df$group)
  missing <- character()
  if (spec$include_all_groups) {
    sources <- all_groups
  } else {
    requested <- spec$include_groups %||% character()
    missing <- setdiff(requested, all_groups)
    sources <- intersect(requested, all_groups)
  }
  sources <- setdiff(sources, spec$exclude_groups %||% character())
  if (length(missing) > 0) {
    if (spec$strict) {
      cli::cli_abort(
        "Aggregate group {.val {spec$group}} references missing group{?s} {.val {missing}}."
      )
    }
    cli::cli_warn(
      "Aggregate group {.val {spec$group}} skips missing group{?s} {.val {missing}}."
    )
  }
  # One source group is a rename, not a total, and it would appear in every
  # table as a duplicate row of the group it copies. Two is the point.
  if (length(sources) < 2) {
    cli::cli_warn(c(
      "Aggregate group {.val {spec$group}} resolved to
       {length(sources)} source group{?s}; skipping.",
      i = "Pooling needs at least two groups -- a single-group fit has
           nothing to total over."
    ))
    return(character())
  }
  sources
}

#' Append synthetic aggregate reporting groups to a draws frame
#'
#' The group-axis counterpart to [add_aggregate_units()]: pools two or more
#' groups into one synthetic group, present for reporting only and never fed
#' back into the fit. Use it to report the implied total across the groups a
#' model was fit on -- the overall effect behind a by-race or by-age split.
#'
#' Each spec aggregates from the *original* frame, so chaining multiple specs
#' never double-counts. An aggregate colliding with an existing group is an
#' error unless the spec sets `overwrite`. A spec resolving to fewer than two
#' source groups is skipped with a warning: pooling one group would just
#' duplicate it under a new name.
#'
#' Applied *after* [add_aggregate_units()] on the same frame, the pooled group
#' also covers the pooled units, giving the "all treated units, all groups"
#' cell. The two orders agree, since both axes pool `mu` by log-sum-exp.
#'
#' The pooled effect is `sum(exp(mu_treated)) / sum(exp(mu))` over the member
#' groups: the total weighted by each group's own exposure, not the mean of
#' the per-group percent changes. Member groups must therefore share a
#' denominator definition for the summed exposure to mean anything -- pooling
#' births-per-population across races is fine, pooling a rate per population
#' with a rate per births is not.
#'
#' @param draws A `countsynth_draws` frame.
#' @param specs A list of [countsynth_aggregate_group()] specs.
#' @return The draws frame with aggregate-group rows appended.
#' @seealso [add_aggregate_units()] to pool the other axis of the cell.
#' @export
add_aggregate_groups <- function(draws, specs) {
  apply_aggregate_specs(
    draws, specs,
    axis = "group", spec_class = "countsynth_aggregate_group",
    resolve_sources = source_groups_for_spec, attr_name = "aggregate_groups"
  )
}

#' Names of the synthetic aggregate groups in a draws frame
#'
#' Reads the `aggregate_groups` attribute stamped by [add_aggregate_groups()],
#' which survives row subsetting and `dplyr` verbs. Returns `character(0)` for
#' a frame that never had any.
#'
#' @param draws A `countsynth_draws` frame.
#' @return Character vector of group names.
#' @export
aggregate_group_names <- function(draws) {
  attr(draws, "aggregate_groups") %||% character()
}
