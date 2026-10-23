# Synthetic aggregate reporting rows. Port of aggregate_units.py: appended to
# the draws frame for reporting/PPC only, never fed back into training.
#
# Two axes of the (unit, time, group) cell can be pooled, and the operation is
# the same either way -- only the key being collapsed differs. Pooling over
# units answers "what happened across all treated states?"; pooling over
# groups (see add_aggregate_groups()) answers "what is the implied total
# across all race groups?". Over whichever axis: outcome/ypred/denominator are
# summed, treatment is max'd, period boundaries carried, and mu/mu_treated
# pooled by log-sum-exp -- the correct pooling for log-count columns, since
# every estimand downstream sums exp(mu).
#
# The two commute (log-sum-exp is associative), so applying both yields a
# well-defined cross cell: the pooled group of the pooled unit.
#
# A pooled unit's max'd treatment marks it exposed from its first member's
# exposure on, which is right for drawing it but wrong for its effect when
# members adopt at different times. Effects are therefore computed from the
# members' own exposed cells instead (post_treatment_rows()), using the member
# sets recorded in the `aggregate_unit_sources` attribute.

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

# Grouped reductions here have to stay in C. Pooling the group axis keys on
# (draw, time, unit), so `include_all_groups` over a 10k-draw, 141-unit,
# 20-period panel asks for 28.6M groups of two rows each: a grouped
# summarise then calls an R-level log-sum-exp 57M times and carries a group
# index costing ~92 bytes per group -- 2.6 GB before a single summary is
# computed, on top of a full copy of the frame. That is what turned a report
# into a 3.3-hour OOM. rowsum() and order() do the same reductions in one C
# pass each.

# max() per group, with max()'s own NA semantics: a group holding any NA
# gives NA. `gid` is never NA, so ordering on (gid, x) sorts NAs last *within*
# each group's block, leaving each group's maximum -- or its NA -- last.
group_max <- function(x, gid, ng) {
  o <- order(gid, x, method = "radix", na.last = TRUE)
  x[o][cumsum(tabulate(gid, ng))]
}

# rowsum() answers with a one-column matrix whose rownames are the group ids;
# unname() keeps those ids from riding along as names on every pooled column.
group_sum <- function(x, gid, ng) {
  unname(rowsum(x, gid, reorder = TRUE)[, 1L])
}

# logsumexp() per group, NAs dropped and an all-NA group giving NA. Shifting
# by each group's own maximum keeps the exponentials in range exactly as
# logsumexp() does, and a non-finite maximum passes straight through the same
# way (all -Inf gives -Inf, any +Inf gives +Inf).
group_logsumexp <- function(x, gid, ng) {
  na <- is.na(x)
  finite_x <- x
  finite_x[na] <- -Inf
  m <- group_max(finite_x, gid, ng)
  shift <- m
  shift[!is.finite(shift)] <- 0
  e <- exp(x - shift[gid])
  e[na] <- 0
  total <- group_sum(e, gid, ng)
  n_obs <- group_sum(as.numeric(!na), gid, ng)
  out <- log(total) + shift
  out[n_obs == 0] <- NA_real_
  out
}

# A single source row's suppressed/missing count would otherwise NA out the
# whole summed cell (and everything downstream that sums it), even when every
# other source has real data -- e.g. one small state missing a subgroup blanks
# the entire "all treated states" total for every period it's missing. Fill
# it with the model's factual expectation instead: mu_treated already equals
# mu on control cells and mu + te on exposed ones, so exp(mu_treated) is
# E[Y(1)] where the cell is exposed and E[Y(0)] otherwise -- whichever
# actually happened there -- and flag the cell in `outcome_imputed` so
# reporting can mark it. The fill must be constant across draws (the
# posterior mean, not the per-draw mu_treated): outcome is supposed to be
# fixed data, and downstream code (compute_quantiles) groups by it to carry it
# through summaries, assuming it never varies within a (unit, time, group)
# cell -- a per-draw fill would silently split every draw into its own group.
#
# The flag is OR'd against any already on the frame rather than overwriting
# it: pooling the second axis sees a cell the first pass already filled,
# whose outcome is no longer NA, and the imputation has to stay visible.
impute_missing_outcomes <- function(sub) {
  missing <- FALSE
  if ("outcome" %in% names(sub)) {
    missing <- is.na(sub$outcome)
    if ("mu_treated" %in% names(sub) && any(missing)) {
      cell_key <- intersect(c("unit", "time", "group"), names(sub))
      cid <- vctrs::vec_group_id(sub[cell_key])
      n_cells <- attr(cid, "n")
      fill <- group_sum(exp(sub$mu_treated), cid, n_cells) /
        tabulate(cid, n_cells)
      sub$outcome[missing] <- fill[cid[missing]]
    }
  }
  prior <- if ("outcome_imputed" %in% names(sub)) {
    !is.na(sub$outcome_imputed) & sub$outcome_imputed
  } else {
    FALSE
  }
  sub$outcome_imputed <- rep_len(missing | prior, nrow(sub))
  sub
}

# Collapse `axis` ("unit" or "group") over `sources` into one synthetic row
# per remaining cell key. Every other key -- draw, time, and the opposite
# axis -- is preserved, so a spec on one axis never disturbs the other.
aggregate_over <- function(source_df, axis, new_name, sources) {
  keep <- source_df[[axis]] %in% sources
  if (!any(keep)) {
    return(NULL)
  }
  kept_axis <- setdiff(c("unit", "group"), axis)
  group_cols <- intersect(
    c(".draw", ".chain", ".iteration", "time", kept_axis), names(source_df)
  )
  # Only the columns the reduction actually reads. Narrowing matters at
  # report scale: when every level of the axis is a source (the
  # `include_all_*` specs), `keep` is all TRUE and this becomes a column
  # selection that copies no rows at all, where subsetting the frame
  # duplicated it outright.
  value_cols <- intersect(
    c("outcome", "ypred", "denominator", "treatment", "mu", "mu_treated",
      "outcome_imputed", "start_date", "end_date"),
    names(source_df)
  )
  needed <- unique(c(group_cols, axis, value_cols))
  sub <- if (all(keep)) source_df[needed] else source_df[keep, needed]
  has <- function(col) col %in% names(sub)
  sub <- impute_missing_outcomes(sub)
  imputed <- sub$outcome_imputed

  gid <- vctrs::vec_group_id(sub[group_cols])
  ng <- attr(gid, "n")
  # vec_group_id() numbers groups by first appearance and rowsum(reorder =
  # TRUE) returns them in increasing id order, so every column below lands in
  # the same order as these keys -- the frame's own row order, rather than
  # the sorted order a grouped summarise would have imposed.
  first_row <- match(seq_len(ng), gid)
  agg <- sub[first_row, group_cols]
  for (col in intersect(c("outcome", "ypred", "denominator"), names(sub))) {
    agg[[col]] <- group_sum(sub[[col]], gid, ng)
  }
  if (has("treatment")) {
    agg$treatment <- group_max(sub$treatment, gid, ng)
  }
  agg$outcome_imputed <- group_sum(as.numeric(imputed), gid, ng) > 0
  for (col in intersect(c("start_date", "end_date"), names(sub))) {
    agg[[col]] <- sub[[col]][first_row]
  }
  for (col in intersect(c("mu", "mu_treated"), names(sub))) {
    agg[[col]] <- group_logsumexp(sub[[col]], gid, ng)
  }
  agg[[axis]] <- new_name
  ordered_cols <- intersect(names(source_df), names(agg))
  extra_cols <- setdiff(names(agg), ordered_cols)
  agg[c(ordered_cols, extra_cols)]
}

# Shared driver behind add_aggregate_units() / add_aggregate_groups(). Every
# spec aggregates from the *original* frame, so chaining specs on one axis
# never double-counts; `resolve_sources` is the axis-specific selector.
# `sources_attr`, when given, names an attribute recording which members each
# aggregate pooled (see post_treatment_rows() for why units need it).
apply_aggregate_specs <- function(draws, specs, axis, spec_class,
                                  resolve_sources, attr_name,
                                  sources_attr = NULL) {
  if (is.null(specs) || length(specs) == 0) {
    return(draws)
  }
  checkmate::assert_list(specs, types = spec_class)
  # aggregate_units/aggregate_groups are carried through too: bind_rows drops
  # attributes, and pooling one axis must not forget what the other stamped.
  attr_names <- c(
    "groups", "units", "times", "has_denominator",
    "aggregate_units", "aggregate_groups", "aggregate_unit_sources"
  )
  attrs <- stats::setNames(lapply(attr_names, function(a) attr(draws, a)), attr_names)
  attrs <- attrs[!vapply(attrs, is.null, logical(1))]
  source_df <- draws
  result_df <- draws
  existing <- unique(draws[[axis]])
  aggregate_frames <- list()
  source_sets <- list()
  replaced <- character()

  for (spec in specs) {
    name <- spec[[axis]]
    if (name %in% existing && !spec$overwrite) {
      cli::cli_abort(
        "Aggregate {axis} {.val {name}} collides with an existing {axis};
         set {.field overwrite} = TRUE to allow replacement."
      )
    }
    if (spec$overwrite) {
      result_df <- result_df[result_df[[axis]] != name, ]
      aggregate_frames[[name]] <- NULL
      source_sets[[name]] <- NULL
      replaced <- c(replaced, name)
      existing <- setdiff(existing, name)
    }
    sources <- resolve_sources(spec, source_df)
    if (length(sources) == 0) {
      next
    }
    aggregate <- aggregate_over(source_df, axis, name, sources)
    if (is.null(aggregate)) {
      next
    }
    aggregate_frames[[name]] <- aggregate
    source_sets[[name]] <- sources
    existing <- c(existing, name)
  }

  if (length(aggregate_frames) > 0 && !"outcome_imputed" %in% names(result_df)) {
    result_df$outcome_imputed <- FALSE
  }
  out <- dplyr::bind_rows(c(list(result_df), unname(aggregate_frames)))
  for (nm in names(attrs)) {
    attr(out, nm) <- attrs[[nm]]
  }
  attr(out, attr_name) <- unique(c(
    attrs[[attr_name]] %||% character(), names(aggregate_frames)
  ))
  if (!is.null(sources_attr)) {
    all_sources <- attrs[[sources_attr]] %||% list()
    all_sources[intersect(replaced, names(all_sources))] <- NULL
    all_sources[names(source_sets)] <- source_sets
    attr(out, sources_attr) <- all_sources
  }
  if (!inherits(out, "countsynth_draws")) {
    class(out) <- c("countsynth_draws", class(out))
  }
  out
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
#' A pooled unit's rows total every member at every time, and its `treatment`
#' is 1 wherever any member's is -- the calendar-time series the fit and gap
#' figures draw. Its *effect* respects staggered adoption instead:
#' [countsynth_summary_table()], [countsynth_post_treatment_summary()] and
#' [countsynth_interval_plot()] sum only over the cells where each member is
#' itself exposed, so a member still waiting on its own start date adds
#' nothing to the counterfactual or the exposure. The member sets behind this
#' are kept in an `aggregate_unit_sources` attribute, and the members' own
#' rows must stay in the frame: with one missing, the effect falls back to the
#' pooled unit's own rows, with a warning.
#'
#' @param draws A `countsynth_draws` frame.
#' @param specs A list of [countsynth_aggregate_unit()] specs.
#' @return The draws frame with aggregate-unit rows appended.
#' @seealso [add_aggregate_groups()] to pool the other axis of the cell.
#' @export
add_aggregate_units <- function(draws, specs) {
  apply_aggregate_specs(
    draws, specs,
    axis = "unit", spec_class = "countsynth_aggregate_unit",
    resolve_sources = source_units_for_spec, attr_name = "aggregate_units",
    sources_attr = "aggregate_unit_sources"
  )
}

# Post-treatment rows as the effect estimands read them, for `units` (NULL =
# every unit).
#
# A pooled unit's own rows can't answer this. Its `treatment` is max'd over
# its members, so it counts as exposed from the first member's exposure on,
# and summing its post-treatment rows folds every other member's
# not-yet-exposed periods into the counterfactual and the exposure. The count
# difference survives that (mu_treated == mu on those cells, so they cancel),
# but every ratio is pulled toward zero: with one state exposed ten months
# before the rest, a third of the pooled post-period can be unexposed. So a
# pooled unit's post-treatment rows are rebuilt here from its members' own
# exposed cells, relabeled to the pooled unit -- each member enters on its own
# adoption date, the usual overall effect under staggered adoption. The
# pooled unit's own rows are left alone for the fit and gap figures, where
# the calendar-time total over every member is the series worth drawing.
#
# `source_unit` names the unit each row came from (the unit itself on real
# rows), since relabeled member rows share their (unit, time, group) key.
post_treatment_rows <- function(draws, units = NULL) {
  treated <- !is.na(draws$treatment) & draws$treatment == 1
  wanted <- if (is.null(units)) treated else treated & draws$unit %in% units
  sources <- attr(draws, "aggregate_unit_sources") %||% list()
  pooled <- names(sources)
  if (!is.null(units)) {
    pooled <- intersect(pooled, units)
  }
  rebuilt <- list()
  if (length(pooled) > 0) {
    real_units <- setdiff(unique(draws$unit), names(sources))
    for (u in pooled) {
      members <- sources[[u]]
      absent <- setdiff(members, real_units)
      # A member can be missing because the frame was subset to the pooled
      # unit, or because the spec overwrote a real unit it also pooled. Its
      # rows are gone either way, so there is nothing to rebuild from.
      if (length(absent) > 0) {
        cli::cli_warn(c(
          "Pooled unit {.val {u}} is missing member unit{?s} {.val {absent}}.",
          i = "Its effect falls back to its own rows, which count every member
               as exposed from the first member's exposure on."
        ))
        next
      }
      rows <- impute_missing_outcomes(draws[treated & draws$unit %in% members, ])
      rows$source_unit <- rows$unit
      rows$unit <- rep(u, nrow(rows))
      rebuilt[[u]] <- rows
    }
  }
  post <- draws[wanted & !draws$unit %in% names(rebuilt), ]
  post$source_unit <- post$unit
  if (length(rebuilt) == 0) {
    return(post)
  }
  if (!"outcome_imputed" %in% names(post)) {
    post$outcome_imputed <- rep(FALSE, nrow(post))
  }
  dplyr::bind_rows(c(list(post), unname(rebuilt)))
}

#' Names of the synthetic aggregate units in a draws frame
#'
#' Reads the `aggregate_units` attribute stamped by [add_aggregate_units()],
#' which survives row subsetting and `dplyr` verbs. Returns `character(0)` for
#' a frame that never had any.
#'
#' @param draws A `countsynth_draws` frame.
#' @return Character vector of unit names.
#' @export
aggregate_unit_names <- function(draws) {
  attr(draws, "aggregate_units") %||% character()
}
