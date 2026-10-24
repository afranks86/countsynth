# Data ingestion: CSV -> standardized long panel -> dense model arrays.
# Faithful port of the Python data.py pipeline; after loading, all operations
# use fixed column names: unit, time, group, outcome, denominator, treatment.

DATE_FORMATS_AUTO <- c("%Y-%m-%d", "%m/%d/%y", "%m/%d/%Y", "%d-%m-%Y")

# Stand-in for a validated (outcome = 0, denominator = 0) cell's scaled
# denominator. These cells are excluded from the likelihood entirely (see
# build_model_arrays()), so this value only keeps log(denominator) finite for
# generated-quantities code paths; it is never fit against data.
ZERO_DENOM_PLACEHOLDER <- 1e-8

parse_time_column <- function(x, date_format, time_col) {
  if (inherits(x, "Date")) {
    return(x)
  }
  if (is.numeric(x)) {
    cli::cli_abort(
      "Could not parse time column {.field {time_col}}: numeric values are not
       supported; provide dates."
    )
  }
  formats <- if (identical(date_format, "auto")) {
    DATE_FORMATS_AUTO
  } else {
    date_format
  }
  for (fmt in formats) {
    parsed <- as.Date(x, format = fmt)
    if (!anyNA(parsed[!is.na(x)])) {
      return(parsed)
    }
  }
  cli::cli_abort("Could not parse time column {.field {time_col}}.")
}

resolve_outcomes_from_prefixes <- function(prefixes, available_columns) {
  outcome_prefix <- prefixes$outcome_prefix
  denominator_prefix <- prefixes$denominator_prefix
  include <- prefixes$include

  outcome_cols <- available_columns[startsWith(available_columns, outcome_prefix)]
  if (length(outcome_cols) == 0) {
    cli::cli_abort(
      "No columns found with outcome prefix {.val {outcome_prefix}}."
    )
  }
  labels <- substring(outcome_cols, nchar(outcome_prefix) + 1L)
  keep <- nzchar(labels)
  outcome_cols <- outcome_cols[keep]
  labels <- labels[keep]

  if (!is.null(include)) {
    missing <- setdiff(include, labels)
    if (length(missing) > 0) {
      cli::cli_abort(
        "Missing outcome column{?s} for include label{?s} {.val {missing}}."
      )
    }
    keep <- labels %in% include
    outcome_cols <- outcome_cols[keep]
    labels <- labels[keep]
  }

  outcomes <- lapply(seq_along(outcome_cols), function(i) {
    denom_col <- NULL
    if (!is.null(denominator_prefix)) {
      denom_col <- paste0(denominator_prefix, labels[[i]])
      if (!denom_col %in% available_columns) {
        cli::cli_abort(
          "Missing denominator column for label {.val {labels[[i]]}}
           (expected {.field {denom_col}})."
        )
      }
    }
    list(
      outcome_col = outcome_cols[[i]], denominator_col = denom_col,
      label = labels[[i]]
    )
  })
  outcomes
}

resolve_outcomes <- function(schema, available_columns) {
  if (!is.null(schema$outcomes)) {
    outcomes <- lapply(schema$outcomes, function(o) {
      list(
        outcome_col = o$outcome_col, denominator_col = o$denominator_col,
        label = o$label
      )
    })
  } else {
    outcomes <- resolve_outcomes_from_prefixes(
      schema$outcomes_from_prefixes, available_columns
    )
  }
  labels <- vapply(outcomes, `[[`, character(1), "label")
  if (anyDuplicated(labels)) {
    cli::cli_abort("Duplicate outcome labels: {.val {sort(labels)}}.")
  }
  outcomes
}

validate_and_resolve_total <- function(groups, outcomes, type_spec) {
  defined_labels <- vapply(outcomes, `[[`, character(1), "label")
  if (!"total" %in% groups || "total" %in% defined_labels) {
    return(NULL)
  }
  total_from <- type_spec$total_from
  total_all <- isTRUE(type_spec$total_all)
  if (!is.null(total_from)) {
    if (anyDuplicated(total_from)) {
      cli::cli_abort("total_from contains duplicate labels: {.val {total_from}}.")
    }
    undefined <- setdiff(total_from, defined_labels)
    if (length(undefined) > 0) {
      cli::cli_abort(
        "total_from references undefined outcome label{?s} {.val {undefined}}."
      )
    }
    return(total_from)
  }
  if (total_all) {
    return(defined_labels)
  }
  cli::cli_abort(
    "'total' is not defined as an outcome label in the schema; provide
     {.field total_from} or {.field total_all} in the model type to create it."
  )
}

# Source column behind each group's denominator, for error messages. A
# synthetic "total" built from `total_from` has no single column, so it maps
# to NA and is reported by group name instead.
denominator_columns <- function(outcomes) {
  stats::setNames(
    vapply(outcomes, function(o) o$denominator_col %||% NA_character_, character(1)),
    vapply(outcomes, `[[`, character(1), "label")
  )
}

# Checked on the long frame *after* the date window, unit exclusions and group
# selection have been applied. Only rows the model actually sees matter here:
# validating the raw wide frame instead used to reject a perfectly good run
# because of a row outside `end_date`, an excluded unit, or an outcome column
# belonging to a group that was never modeled.
#
# The denominator is an exposure/offset, not a trial count, so an outcome
# exceeding it is not by itself invalid: the two can be on different scales
# (e.g. a denominator reported in thousands). What's fatal is a mismatch no
# rescaling can explain: a negative denominator, a negative outcome, or a
# positive outcome against a literal zero denominator (no amount of exposure
# is zero and positive at once). A denominator of 0 alongside an observed
# outcome of 0 is not invalid data either: zero exposure deterministically
# implies zero count (MCAR, not informative missingness), so that cell is
# excluded from the likelihood entirely rather than fed a literal log(0) (see
# the `zero_denom` handling in build_model_arrays()). A row whose outcome is
# NA (suppressed/missing) keeps the original all-or-nothing rule, since a
# missing count can't confirm the zero-exposure case above.
validate_denominators <- function(df_long, outcomes) {
  if (!"denominator" %in% names(df_long) || nrow(df_long) == 0) {
    return(invisible(NULL))
  }
  denom <- df_long$denominator
  outcome <- df_long$outcome
  observed <- !is.na(outcome)
  bad <- ifelse(
    observed,
    is.na(denom) | denom < 0 | outcome < 0 | (denom == 0 & outcome > 0),
    is.na(denom) | denom <= 0
  )
  if (!any(bad)) {
    return(invisible(NULL))
  }
  rows <- df_long[bad, ]
  cols <- denominator_columns(outcomes)
  sources <- unique(vapply(unique(rows$group), function(g) {
    col <- if (g %in% names(cols)) cols[[g]] else NA_character_
    if (is.na(col)) sprintf("group %s", g) else sprintf("column %s", col)
  }, character(1)))
  times <- sort(unique(rows$time))
  cli::cli_abort(c(
    "Invalid outcome/denominator values in {sources}: {sum(bad)} row{?s}
     affected (unit{?s}: {.val {unique(rows$unit)}}).",
    i = "Period{?s}: {.val {format(times)}}.",
    i = "Each row needs a non-negative outcome and a non-negative denominator
         (a (0, 0) pair is fine and is treated as a deterministic zero, but a
         positive outcome needs a positive denominator); an unobserved row
         still needs a positive denominator. Rows outside the configured
         date window, excluded units, and unmodeled groups are already
         dropped before this check, so every row named here is one the fit
         would have used."
  ))
}

# Every row carries the range its true count is known to lie in alongside
# the count itself: an observed count is its own range, and a suppressed one
# (NA) lies in `suppression_range`. Pooling rows -- a synthetic total here, a
# time block in summarise_blocks() -- sums the bounds, so a pooled cell with
# a suppressed part keeps the range it really has rather than being mistaken
# for a lone suppressed cell.
wide_to_long <- function(df, schema, outcomes, groups, total_from_labels,
                         suppression_range = DEFAULT_SUPPRESSION_RANGE) {
  defined_labels <- vapply(outcomes, `[[`, character(1), "label")
  needs_total <- "total" %in% groups && !"total" %in% defined_labels
  total_labels <- if (needs_total && !is.null(total_from_labels)) {
    total_from_labels
  } else {
    character()
  }

  keep <- vapply(
    outcomes,
    function(o) o$label %in% groups || o$label %in% total_labels,
    logical(1)
  )
  if (!any(keep)) {
    cli::cli_abort("No matching outcomes found for groups: {.val {groups}}.")
  }

  long_parts <- lapply(outcomes[keep], function(o) {
    part <- tibble::tibble(
      unit = df[[schema$unit_col]],
      time = df[[schema$time_col]],
      treatment = df[[schema$treatment_col]],
      group = o$label,
      outcome = as.numeric(df[[o$outcome_col]])
    )
    suppressed <- is.na(part$outcome)
    part$outcome_lower <- ifelse(suppressed, suppression_range[[1]], part$outcome)
    part$outcome_upper <- ifelse(suppressed, suppression_range[[2]], part$outcome)
    if (!is.null(o$denominator_col) && o$denominator_col %in% names(df)) {
      part$denominator <- as.numeric(df[[o$denominator_col]])
    }
    part
  })
  df_long <- dplyr::bind_rows(long_parts)

  if (needs_total && length(total_labels) > 0) {
    df_for_total <- dplyr::filter(df_long, .data$group %in% total_labels)
    has_denom <- "denominator" %in% names(df_for_total)
    df_total <- df_for_total |>
      dplyr::group_by(.data$unit, .data$time, .data$treatment) |>
      dplyr::summarise(
        outcome = sum(.data$outcome),
        outcome_lower = sum(.data$outcome_lower),
        outcome_upper = sum(.data$outcome_upper),
        denominator = if (has_denom) sum(.data$denominator) else 1,
        .groups = "drop"
      ) |>
      dplyr::mutate(group = "total")
    df_long <- dplyr::filter(df_long, .data$group %in% groups)
    df_long <- dplyr::bind_rows(df_long, df_total)
  } else {
    df_long <- dplyr::filter(df_long, .data$group %in% groups)
  }

  dplyr::arrange(df_long, .data$unit, .data$time, .data$group)
}

# Cell-level aggregation rule, shared by both binning strategies: counts sum,
# and treatment is max'd so a partly exposed block counts as exposed. The
# denominator follows what it measures. A stock (population standing through
# the period) averages, and reporting turns it into person-time by weighting
# with the block's length. A flow (births *during* the period) sums, exactly
# as the counts do: averaging it instead left a quarterly deaths-per-birth
# rate at 3x its true value, since the deaths were summed over three months
# against one month's births.
#
# The censoring bounds sum too (see wide_to_long()), so a block pooling an
# observed 50, an observed 48 and one suppressed month comes out NA with the
# range 99..107 -- not the 1..9 a lone suppressed cell would get. `.n_rows`
# counts the input periods each (unit, block, group) actually had, for
# check_block_coverage().
summarise_blocks <- function(df, block_cols, denominator_flow = FALSE) {
  has_denom <- "denominator" %in% names(df)
  has_bounds <- all(c("outcome_lower", "outcome_upper") %in% names(df))
  combine_denom <- if (denominator_flow) sum else mean
  out <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(block_cols))) |>
    dplyr::summarise(
      outcome = sum(.data$outcome),
      outcome_lower = if (has_bounds) sum(.data$outcome_lower) else NA_real_,
      outcome_upper = if (has_bounds) sum(.data$outcome_upper) else NA_real_,
      treatment = max(.data$treatment),
      denominator = if (has_denom) combine_denom(.data$denominator) else NA_real_,
      .n_rows = dplyr::n(),
      .groups = "drop"
    )
  if (!has_bounds) {
    out$outcome_lower <- NULL
    out$outcome_upper <- NULL
  }
  # "No denominator" means every cell behaves as if denominator == 1 (see
  # build_model_arrays()), which depends on the `denominator` column being
  # absent entirely, not present-and-1: build_model_arrays() only skips
  # scaling by `denominator_scale` when the column doesn't exist. Fabricating
  # a literal 1 here would otherwise get divided down to 1/denominator_scale
  # downstream, silently rescaling every rate the moment aggregation is on.
  if (!has_denom) {
    out$denominator <- NULL
  }
  out
}

# A block summed over fewer input periods than the panel has in it would
# undercount with nothing to show it: a quarter one unit reported only two
# months of looks like a one-third drop in its rate, and the unbalanced-panel
# check in build_model_arrays() never sees it, because the block itself
# exists. So an incomplete block is an unbalanced panel in its own right --
# an error unless the config allows one, and otherwise dropped, which leaves
# the cell structurally absent and out of the likelihood. `expected` is the
# number of input periods the panel has in each row's block, so a bin the
# date window cuts short for every unit is not mistaken for one unit's gap.
check_block_coverage <- function(df_agg, expected, allow_unbalanced_panel) {
  short <- df_agg$.n_rows < expected
  df_agg$.n_rows <- NULL
  if (!any(short)) {
    return(list(df = df_agg, note = character()))
  }
  gaps <- unique(df_agg[short, c("unit", "time")])
  examples <- utils::head(sprintf("%s %s", gaps$unit, format(gaps$time)), 3)
  if (!allow_unbalanced_panel) {
    cli::cli_abort(c(
      "Unbalanced panel: {nrow(gaps)} unit-block{?s} {?is/are} missing input
       periods inside a time-aggregation block (e.g. {.val {examples}}).",
      i = "Summed over the periods that are present, such a block would
           undercount with nothing to show it.",
      i = "Fill the gaps, aggregate the data yourself before fitting, or set
           {.field allow_unbalanced_panel} = TRUE to leave these blocks out of
           the likelihood."
    ))
  }
  note <- cli::format_inline(
    "{nrow(gaps)} unit-block{?s} {?was/were} missing input periods and
     {?was/were} left out of the likelihood (e.g. {.val {examples}})."
  )
  list(df = df_agg[!short, ], note = note)
}

# Coverage of each distinct input period: from its own start to the day
# before the next one begins. Only the last period has no successor, and a
# gap in the panel (a year missing for every unit) would otherwise stretch
# the period before it across the gap, so each period is also capped at the
# panel's usual step -- in whole months when every time falls on the first
# of a month, since months differ in length, and in days otherwise. Non-Date
# times (an integer index) get an exclusive right edge instead.
period_bounds <- function(times) {
  times <- sort(unique(times))
  nt <- length(times)
  modal <- function(x) as.numeric(names(which.max(table(x))))
  if (!inherits(times, "Date")) {
    step <- if (nt > 1) modal(diff(times)) else 1
    nxt <- c(times[-1], times[nt] + step)
    return(data.frame(time = times, start = times, end = pmin(nxt, times + step)))
  }
  if (all(format(times, "%d") == "01")) {
    month_index <- as.integer(format(times, "%Y")) * 12L +
      as.integer(format(times, "%m"))
    step <- if (nt > 1) modal(diff(month_index)) else 1
    cap <- add_months(times, step)
  } else {
    step <- if (nt > 1) modal(as.numeric(diff(times))) else 1
    cap <- times + step
  }
  nxt <- c(times[-1], cap[nt])
  data.frame(time = times, start = times, end = pmin(nxt, cap) - 1)
}

# Period bounds for a panel that was not aggregated here, so reporting can
# weight a stock denominator by each period's real length. Without them every
# row counts as one whole year: a bimonthly panel's "person-years" came out
# six times too large and every rate six times too small.
attach_period_bounds <- function(df) {
  if (nrow(df) == 0 || !inherits(df$time, "Date")) {
    return(df)
  }
  bounds <- period_bounds(df$time)
  idx <- match(df$time, bounds$time)
  df$start_date <- bounds$start[idx]
  df$end_date <- bounds$end[idx]
  df
}

# Calendar binning: bins land on calendar boundaries regardless of when the
# panel starts. Assumes monthly-or-finer input, since it bins on year+month.
aggregate_calendar <- function(df, period, denominator_flow = FALSE,
                               allow_unbalanced_panel = FALSE) {
  months_per_period <- switch(period,
    monthly = 1L, bimonthly = 2L, quarterly = 3L, yearly = 12L,
    cli::cli_abort(
      "Unknown aggregation period {.val {period}}; use monthly, bimonthly,
       quarterly, or yearly."
    )
  )
  bin_of <- function(time) {
    month <- as.integer(format(time, "%m"))
    data.frame(
      .year = as.integer(format(time, "%Y")),
      .period = if (period == "yearly") 1L else (month - 1L) %/% months_per_period + 1L
    )
  }
  df[c(".year", ".period")] <- bin_of(df$time)
  df_agg <- summarise_blocks(
    df, c("unit", ".year", ".period", "group"), denominator_flow
  )

  # What each bin actually covers, from the input periods the panel has in
  # it. That is the calendar bin itself for monthly data filling it, but not
  # for a bin the date window cuts short (a quarter holding only October),
  # nor for input coarser than the bin (half-years in quarterly bins), and
  # reporting turns these dates into person-time.
  bounds <- period_bounds(df$time)
  bounds[c(".year", ".period")] <- bin_of(bounds$time)
  bins <- bounds |>
    dplyr::group_by(.data$.year, .data$.period) |>
    dplyr::summarise(
      .expected = dplyr::n(), .start = min(.data$start), .end = max(.data$end),
      .groups = "drop"
    )
  bins$.bin_start <- as.Date(sprintf(
    "%d-%02d-01", bins$.year, (bins$.period - 1L) * months_per_period + 1L
  ))
  bin_end <- add_months(bins$.bin_start, months_per_period) - 1L
  notes <- character()
  overhang <- bins$.end > bin_end
  if (any(overhang)) {
    notes <- c(notes, cli::format_inline(
      "The input periods are longer than a {period} bin, so binning
       combines nothing; each block keeps its input period's dates."
    ))
  }
  partial <- !overhang & (bins$.start > bins$.bin_start | bins$.end < bin_end)
  if (any(partial)) {
    cut <- format(bins$.bin_start[partial])
    notes <- c(notes, cli::format_inline(
      "{sum(partial)} bin{?s} {?is/are} only partly covered by the data
       ({.val {cut}}); {?its/their} shorter exposure is carried in
       {.field start_date}/{.field end_date}."
    ))
  }

  df_agg <- dplyr::left_join(df_agg, bins, by = c(".year", ".period"))
  df_agg$time <- df_agg$.bin_start
  df_agg$start_date <- df_agg$.start
  df_agg$end_date <- df_agg$.end
  coverage <- check_block_coverage(
    df_agg, df_agg$.expected, allow_unbalanced_panel
  )
  df_agg <- coverage$df
  df_agg[c(".year", ".period", ".expected", ".start", ".end", ".bin_start")] <- NULL
  df_agg <- dplyr::arrange(df_agg, .data$unit, .data$time, .data$group)
  attr(df_agg, "aggregation_notes") <- c(notes, coverage$note)
  df_agg
}

# Positional binning: combine every `n` consecutive time points, whatever the
# panel's native resolution is. Blocks are cut from the panel's sorted
# distinct times (not per unit), so every unit lands in the same blocks even
# when some unit is missing a period.
aggregate_n_periods <- function(df, n, denominator_flow = FALSE,
                                allow_unbalanced_panel = FALSE) {
  times <- sort(unique(df$time))
  nt <- length(times)
  n_blocks <- ceiling(nt / n)
  remainder <- nt %% n
  notes <- character()
  if (remainder != 0) {
    notes <- c(notes, cli::format_inline(
      "{nt} time point{?s} {?does/do} not divide evenly into blocks of {n};
       the final block holds {remainder} period{?s}. Its shorter exposure is
       carried in {.field start_date}/{.field end_date}."
    ))
  }
  df$.block <- (match(df$time, times) - 1L) %/% n + 1L
  df_agg <- summarise_blocks(df, c("unit", ".block", "group"), denominator_flow)

  # A block runs from its first input period's start to its last one's end
  # (see period_bounds()).
  bounds <- period_bounds(times)
  first_idx <- (seq_len(n_blocks) - 1L) * n + 1L
  last_idx <- pmin(seq_len(n_blocks) * n, nt)
  blocks <- data.frame(
    .block = seq_len(n_blocks),
    .start = times[first_idx],
    .end = bounds$end[last_idx],
    .expected = last_idx - first_idx + 1L
  )
  df_agg <- dplyr::left_join(df_agg, blocks, by = ".block")
  df_agg$time <- df_agg$.start
  df_agg$start_date <- df_agg$.start
  df_agg$end_date <- df_agg$.end
  coverage <- check_block_coverage(
    df_agg, df_agg$.expected, allow_unbalanced_panel
  )
  df_agg <- coverage$df
  df_agg[c(".block", ".start", ".end", ".expected")] <- NULL
  df_agg <- dplyr::arrange(df_agg, .data$unit, .data$time, .data$group)
  attr(df_agg, "aggregation_notes") <- c(notes, coverage$note)
  df_agg
}

# Time aggregation inside the package is a convenience with sharp edges, so
# every run that uses it says what it did to this data, in one warning,
# rather than leaving those choices to be discovered in the numbers.
warn_time_aggregation <- function(df, spec, denominator_flow, has_denominator,
                                  adjust_for_missingness,
                                  suppression_range = DEFAULT_SUPPRESSION_RANGE) {
  bins <- if (!is.null(spec$n_periods)) {
    sprintf("blocks of %d periods", spec$n_periods)
  } else {
    sprintf("%s bins", spec$period)
  }
  bullets <- c(
    i = "Counts are summed and the treatment indicator is max'd, so a block
         is exposed if any period in it is."
  )
  if (has_denominator) {
    bullets <- c(bullets, i = if (denominator_flow) {
      "Denominators are summed within each block: {.field
       output.denominator_time_unit} is {.val none}, which marks a
       per-period flow such as births."
    } else {
      "Denominators are averaged within each block, which is right for a
       level such as population. If yours is a per-period flow such as
       births, set {.field output.denominator_time_unit} to {.val none} so
       it is summed like the counts."
    })
  }
  pooled <- if (all(c("outcome_lower", "outcome_upper") %in% names(df))) {
    sum(is.na(df$outcome) & (df$outcome_lower != suppression_range[[1]] |
      df$outcome_upper != suppression_range[[2]]))
  } else {
    0L
  }
  if (pooled > 0) {
    bullets <- c(bullets, i = if (adjust_for_missingness) {
      cli::format_inline(
        "{pooled} block{?s} pool{?s/} a suppressed period with others;
         {?it is/each is} censored on {?its/their} summed range rather than
         on {suppression_range[[1]]}..{suppression_range[[2]]}."
      )
    } else {
      cli::format_inline(
        "{pooled} block{?s} pool{?s/} a suppressed period with others and,
         with {.field adjust_for_missingness} off, {?is/are} left out of the
         likelihood whole -- observed periods included."
      )
    })
  }
  notes <- attr(df, "aggregation_notes") %||% character()
  # The notes are already formatted; keep cli from reading a brace in a unit
  # name as markup.
  notes <- gsub("}", "}}", gsub("{", "{{", notes, fixed = TRUE), fixed = TRUE)
  bullets <- c(bullets, stats::setNames(notes, rep("i", length(notes))))
  cli::cli_warn(c(
    "!" = "Time aggregation is on ({bins}). Aggregating inside countsynth is
           easy to get subtly wrong; where you can, aggregate the data before
           fitting and leave {.field time_aggregation} off. See
           {.help countsynth::countsynth_time_aggregation} for what to check.",
    bullets
  ))
}

aggregate_temporal <- function(df, spec, denominator_flow = FALSE,
                               allow_unbalanced_panel = FALSE) {
  if (!is.null(spec$n_periods)) {
    aggregate_n_periods(
      df, spec$n_periods, denominator_flow, allow_unbalanced_panel
    )
  } else {
    aggregate_calendar(
      df, spec$period, denominator_flow, allow_unbalanced_panel
    )
  }
}

# Add whole months to a Date (day-of-month is always 1 here).
add_months <- function(dates, months) {
  lt <- as.POSIXlt(dates)
  lt$mon <- lt$mon + months
  as.Date(lt)
}

#' Load and prepare panel data for modeling
#'
#' Reads the config's input CSV, resolves outcomes, reshapes to a standardized
#' long panel (`unit`, `time`, `group`, `outcome`, `denominator`,
#' `treatment`), applies date filtering / unit exclusion / temporal
#' aggregation, and builds the dense `(K, D, N)` model arrays.
#'
#' @param config A [countsynth_config()] object.
#' @param type Name of the model type in `config$model$types` to prepare data
#'   for. Defaults to the first type.
#' @param df Optional data frame to use instead of reading
#'   `config$input_file`.
#' @return A `countsynth_data` object: list with `Y`, `denominators`,
#'   `control_idx_array`, `missing_idx_array`, `excluded_idx_array` (all
#'   `(K, D, N)` arrays), `groups`, `units`, `times`, `type`, and the
#'   standardized long frame `df`.
#' @export
countsynth_data <- function(config, type = NULL, df = NULL) {
  checkmate::assert_class(config, "countsynth_config")
  if (length(config$model$types) == 0) {
    cli::cli_abort("Config must define at least one entry in {.field model.types}.")
  }
  type <- type %||% names(config$model$types)[[1]]
  checkmate::assert_choice(type, names(config$model$types))
  type_spec <- config$model$types[[type]]
  groups <- type_spec$groups

  if (is.null(df)) {
    checkmate::assert_file_exists(config$input_file)
    df <- utils::read.csv(config$input_file, check.names = FALSE)
  }
  schema <- config$schema
  outcomes <- resolve_outcomes(schema, names(df))

  required <- unique(c(
    schema$unit_col, schema$time_col, schema$treatment_col,
    vapply(outcomes, `[[`, character(1), "outcome_col"),
    unlist(lapply(outcomes, `[[`, "denominator_col"))
  ))
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    cli::cli_abort("Missing column{?s} in input data: {.val {sort(missing_cols)}}.")
  }
  df[[schema$time_col]] <- parse_time_column(
    df[[schema$time_col]], config$date_format, schema$time_col
  )

  total_from_labels <- validate_and_resolve_total(groups, outcomes, type_spec)
  suppression_range <- config$model$suppression_range %||% DEFAULT_SUPPRESSION_RANGE
  df_long <- wide_to_long(
    df, schema, outcomes, groups, total_from_labels, suppression_range
  )

  dupes <- duplicated(df_long[c("group", "unit", "time")])
  if (any(dupes)) {
    sample_units <- utils::head(sort(unique(df_long$unit[dupes])), 3)
    cli::cli_abort(
      "Duplicate group/unit/time combinations: {sum(dupes)} row{?s} affected
       (e.g. units: {.val {sample_units}})."
    )
  }

  if (!is.null(config$start_date)) {
    df_long <- df_long[df_long$time >= as.Date(config$start_date), ]
  }
  if (!is.null(config$end_date)) {
    df_long <- df_long[df_long$time < as.Date(config$end_date), ]
  }

  if (!is.null(type_spec$exclude_units) && length(type_spec$exclude_units) > 0) {
    df_long <- df_long[!df_long$unit %in% type_spec$exclude_units, ]
  }

  validate_denominators(df_long, outcomes)

  if (config$time_aggregation$enabled) {
    # `denominator_time_unit: none` is how a config says its denominator is a
    # per-period flow rather than a standing level (see
    # countsynth_output_opts()), and that one setting decides both how the
    # denominator is weighted in reporting and how it combines here.
    denominator_flow <- identical(config$output$denominator_time_unit, "none")
    df_long <- aggregate_temporal(
      df_long, config$time_aggregation,
      denominator_flow = denominator_flow,
      allow_unbalanced_panel = config$allow_unbalanced_panel
    )
    warn_time_aggregation(
      df_long, config$time_aggregation, denominator_flow,
      has_denominator = "denominator" %in% names(df_long),
      adjust_for_missingness = config$model$adjust_for_missingness,
      suppression_range = suppression_range
    )
  } else {
    df_long <- attach_period_bounds(df_long)
  }

  arrays <- build_model_arrays(
    df_long, groups,
    allow_unbalanced_panel = config$allow_unbalanced_panel,
    suppression_range = suppression_range
  )
  arrays$df <- df_long
  arrays$type <- type
  new_countsynth_class(arrays, "countsynth_data")
}

#' Convert a standardized long frame into dense (K, D, N) model arrays
#'
#' @param df Long data frame with columns `unit`, `time`, `group`, `outcome`,
#'   optionally `denominator`, and `treatment`.
#' @param groups Group labels in the desired order (defines the K dimension).
#' @param denominator_scale Denominators are divided by this (default `1e4`,
#'   i.e. rates per 10k).
#' @param allow_unbalanced_panel If `FALSE`, error when any (group, unit,
#'   time) cell has no row; if `TRUE`, leave such cells out of the likelihood.
#'   A cell with no row was never reported, which is not the same as a
#'   suppressed count, so it gets no censoring term either.
#' @param suppression_range The range a suppressed (`NA`) outcome lies in,
#'   for rows without their own `outcome_lower`/`outcome_upper` bounds.
#' @return List with `Y`, `denominators`, `control_idx_array`,
#'   `missing_idx_array`, `excluded_idx_array`, `cens_lower`, `cens_upper`
#'   (each suppressed cell's range, `NA` elsewhere), `suppression_range`,
#'   `groups`, `units`, `times`.
#' @export
build_model_arrays <- function(df, groups, denominator_scale = 1e4,
                               allow_unbalanced_panel = FALSE,
                               suppression_range = DEFAULT_SUPPRESSION_RANGE) {
  df <- df[df$group %in% groups, ]
  df <- dplyr::arrange(df, .data$unit, .data$time, .data$group)

  units <- sort(unique(df$unit))
  times <- sort(unique(df$time))
  K <- length(groups)
  D <- length(units)
  N <- length(times)

  Y <- array(0, dim = c(K, D, N))
  denominators <- array(1, dim = c(K, D, N))
  control_idx <- array(TRUE, dim = c(K, D, N))
  missing_idx <- array(FALSE, dim = c(K, D, N))
  excluded_idx <- array(FALSE, dim = c(K, D, N))
  filled <- array(FALSE, dim = c(K, D, N))

  k_idx <- match(df$group, groups)
  d_idx <- match(df$unit, units)
  n_idx <- match(df$time, times)
  cells <- cbind(k_idx, d_idx, n_idx)

  outcome_missing <- is.na(df$outcome)
  Y[cells[!outcome_missing, , drop = FALSE]] <-
    as.numeric(df$outcome[!outcome_missing])
  missing_idx[cells[outcome_missing, , drop = FALSE]] <- TRUE

  # Each suppressed cell's range: suppression_range for a cell suppressed on
  # its own, wider for one pooling suppressed and observed counts (see
  # wide_to_long()). A hand-built frame without bounds is all lone cells.
  has_bounds <- all(c("outcome_lower", "outcome_upper") %in% names(df))
  lower <- if (has_bounds) df$outcome_lower else rep(suppression_range[[1]], nrow(df))
  upper <- if (has_bounds) df$outcome_upper else rep(suppression_range[[2]], nrow(df))
  cens_lower <- array(NA_real_, dim = c(K, D, N))
  cens_upper <- array(NA_real_, dim = c(K, D, N))
  cens_lower[cells[outcome_missing, , drop = FALSE]] <- lower[outcome_missing]
  cens_upper[cells[outcome_missing, , drop = FALSE]] <- upper[outcome_missing]

  if ("denominator" %in% names(df)) {
    denom <- as.numeric(df$denominator)
    denom_ok <- !is.na(denom) & denom > 0
    denominators[cells[denom_ok, , drop = FALSE]] <-
      denom[denom_ok] / denominator_scale

    # Zero exposure deterministically implies zero count (validate_denominators()
    # only lets a zero denominator through when the observed outcome is also
    # zero): exclude these cells from the likelihood entirely rather than
    # feed Stan a literal log(0). The tiny placeholder keeps mu_ctrl finite
    # for the counterfactual posterior predictive (ypred); it never reaches
    # the likelihood since these cells are dropped from
    # obs_cell/cens_cell/icens_cell in stan_data_joint()/stan_data_stage2().
    zero_denom <- !is.na(denom) & denom == 0
    denominators[cells[zero_denom, , drop = FALSE]] <- ZERO_DENOM_PLACEHOLDER
    excluded_idx[cells[zero_denom, , drop = FALSE]] <- TRUE
  }

  control_idx[cells] <- df$treatment == 0
  filled[cells] <- TRUE

  n_absent <- sum(!filled)
  if (n_absent > 0) {
    if (!allow_unbalanced_panel) {
      cli::cli_abort(
        "Unbalanced panel: {n_absent} of {K * D * N} cells are missing.
         Set {.field allow_unbalanced_panel} = TRUE to allow."
      )
    }
    # Absent, not suppressed: a cell nobody reported carries no information
    # about its count, so it is excluded outright. Marking it missing gave it
    # a "count was 1..9" censoring term, which drags the fit toward a small
    # count wherever a unit simply has no data.
    excluded_idx <- excluded_idx | !filled
    cli::cli_warn(
      "Unbalanced panel: {n_absent} of {K * D * N} cells are structurally
       absent (left out of the likelihood)."
    )
  }

  list(
    Y = Y, denominators = denominators, control_idx_array = control_idx,
    missing_idx_array = missing_idx, excluded_idx_array = excluded_idx,
    cens_lower = cens_lower, cens_upper = cens_upper,
    suppression_range = as.integer(suppression_range),
    groups = groups, units = units, times = times
  )
}

#' @export
print.countsynth_data <- function(x, ...) {
  K <- length(x$groups)
  D <- length(x$units)
  N <- length(x$times)
  n_exposed <- sum(!x$control_idx_array)
  cli::cli_h1("countsynth data ({x$type})")
  cli::cli_li("{K} group{?s}: {.val {x$groups}}")
  cli::cli_li("{D} unit{?s}, {N} time period{?s}")
  cli::cli_li("{n_exposed} exposed cell{?s}, {sum(x$missing_idx_array)} suppressed cell{?s}")
  if (!is.null(x$excluded_idx_array) && any(x$excluded_idx_array)) {
    cli::cli_li(
      "{sum(x$excluded_idx_array)} cell{?s} left out of the likelihood (zero
       exposure, or absent from the panel)"
    )
  }
  invisible(x)
}
