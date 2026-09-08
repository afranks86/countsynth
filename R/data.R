# Data ingestion: CSV -> standardized long panel -> dense model arrays.
# Faithful port of the Python data.py pipeline; after loading, all operations
# use fixed column names: unit, time, group, outcome, denominator, treatment.

DATE_FORMATS_AUTO <- c("%Y-%m-%d", "%m/%d/%y", "%m/%d/%Y", "%d-%m-%Y")

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
# selection have been applied. A zero or missing denominator only matters for
# rows the model actually sees: stan_data computes log(denominator), so a zero
# inside the window is genuinely fatal (-Inf), while one in a row the config
# already discarded is not a problem at all. Validating the raw wide frame
# instead used to reject a perfectly good run because of a row outside
# `end_date`, an excluded unit, or an outcome column belonging to a group that
# was never modeled.
validate_denominators <- function(df_long, outcomes) {
  if (!"denominator" %in% names(df_long) || nrow(df_long) == 0) {
    return(invisible(NULL))
  }
  bad <- is.na(df_long$denominator) | df_long$denominator <= 0
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
    "Missing or non-positive denominator in {sources}: {sum(bad)} row{?s}
     affected (unit{?s}: {.val {unique(rows$unit)}}).",
    i = "Period{?s}: {.val {format(times)}}.",
    i = "The model uses log(denominator), so these cells have no defined
         exposure. Rows outside the configured date window, excluded units,
         and unmodeled groups are already dropped before this check, so every
         row named here is one the fit would have used."
  ))
}

wide_to_long <- function(df, schema, outcomes, groups, total_from_labels) {
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
# denominators average (they are stocks, not flows), and treatment is max'd so
# a partly exposed block counts as exposed.
summarise_blocks <- function(df, block_cols) {
  has_denom <- "denominator" %in% names(df)
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(block_cols))) |>
    dplyr::summarise(
      outcome = sum(.data$outcome),
      treatment = max(.data$treatment),
      denominator = if (has_denom) mean(.data$denominator) else 1,
      .groups = "drop"
    )
}

# Calendar binning: bins land on calendar boundaries regardless of when the
# panel starts. Assumes monthly-or-finer input, since it bins on year+month.
aggregate_calendar <- function(df, period) {
  months_per_period <- switch(period,
    monthly = 1L, bimonthly = 2L, quarterly = 3L, yearly = 12L,
    cli::cli_abort(
      "Unknown aggregation period {.val {period}}; use monthly, bimonthly,
       quarterly, or yearly."
    )
  )
  month <- as.integer(format(df$time, "%m"))
  df$.year <- as.integer(format(df$time, "%Y"))
  df$.period <- if (period == "yearly") {
    rep(1L, nrow(df))
  } else {
    (month - 1L) %/% months_per_period + 1L
  }

  df_agg <- summarise_blocks(df, c("unit", ".year", ".period", "group"))

  first_month <- (df_agg$.period - 1L) * months_per_period + 1L
  df_agg$time <- as.Date(sprintf("%d-%02d-01", df_agg$.year, first_month))
  df_agg$start_date <- df_agg$time
  df_agg$end_date <- add_months(df_agg$time, months_per_period) - 1L

  df_agg$.year <- NULL
  df_agg$.period <- NULL
  dplyr::arrange(df_agg, .data$unit, .data$time, .data$group)
}

# Positional binning: combine every `n` consecutive time points, whatever the
# panel's native resolution is. Blocks are cut from the panel's sorted
# distinct times (not per unit), so every unit lands in the same blocks even
# when some unit is missing a period.
aggregate_n_periods <- function(df, n) {
  times <- sort(unique(df$time))
  n_blocks <- ceiling(length(times) / n)
  remainder <- length(times) %% n
  if (remainder != 0) {
    cli::cli_warn(
      "{length(times)} time point{?s} does not divide evenly into blocks of
       {n}; the final block holds {remainder} period{?s}. Its shorter
       exposure is carried in {.field start_date}/{.field end_date}."
    )
  }
  # A block runs up to the start of the next input period, which is exact and
  # needs no assumption about period width -- important for calendar months,
  # which are not all the same length. Only the very last period has no
  # successor; extend it by the final observed gap.
  nt <- length(times)
  next_start <- c(
    times[-1],
    times[nt] + if (nt > 1) (times[nt] - times[nt - 1]) else 0
  )
  df$.block <- (match(df$time, times) - 1L) %/% n + 1L

  df_agg <- summarise_blocks(df, c("unit", ".block", "group"))
  is_date <- inherits(times, "Date")
  last_idx <- pmin(seq_len(n_blocks) * n, nt)
  bounds <- data.frame(
    .block = seq_len(n_blocks),
    .start = times[(seq_len(n_blocks) - 1L) * n + 1L],
    # Inclusive last day for dates; exclusive right edge otherwise.
    .end = next_start[last_idx] - if (is_date) 1 else 0
  )
  df_agg <- dplyr::left_join(df_agg, bounds, by = ".block")

  df_agg$time <- df_agg$.start
  df_agg$start_date <- df_agg$.start
  df_agg$end_date <- df_agg$.end

  df_agg$.block <- NULL
  df_agg$.start <- NULL
  df_agg$.end <- NULL
  dplyr::arrange(df_agg, .data$unit, .data$time, .data$group)
}

aggregate_temporal <- function(df, spec) {
  if (!is.null(spec$n_periods)) {
    aggregate_n_periods(df, spec$n_periods)
  } else {
    aggregate_calendar(df, spec$period)
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
#' @param config A [bpnmf_config()] object.
#' @param type Name of the model type in `config$model$types` to prepare data
#'   for. Defaults to the first type.
#' @param df Optional data frame to use instead of reading
#'   `config$input_file`.
#' @return A `bpnmf_data` object: list with `Y`, `denominators`,
#'   `control_idx_array`, `missing_idx_array` (all `(K, D, N)` arrays),
#'   `groups`, `units`, `times`, `type`, and the standardized long frame `df`.
#' @export
bpnmf_data <- function(config, type = NULL, df = NULL) {
  checkmate::assert_class(config, "bpnmf_config")
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
  df_long <- wide_to_long(df, schema, outcomes, groups, total_from_labels)

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
    df_long <- aggregate_temporal(df_long, config$time_aggregation)
  }

  arrays <- build_model_arrays(
    df_long, groups,
    allow_unbalanced_panel = config$allow_unbalanced_panel
  )
  arrays$df <- df_long
  arrays$type <- type
  new_bpnmf_class(arrays, "bpnmf_data")
}

#' Convert a standardized long frame into dense (K, D, N) model arrays
#'
#' @param df Long data frame with columns `unit`, `time`, `group`, `outcome`,
#'   optionally `denominator`, and `treatment`.
#' @param groups Group labels in the desired order (defines the K dimension).
#' @param denominator_scale Denominators are divided by this (default `1e4`,
#'   i.e. rates per 10k).
#' @param allow_unbalanced_panel If `FALSE`, error when any (group, unit,
#'   time) cell has no row; if `TRUE`, mark such cells missing.
#' @return List with `Y`, `denominators`, `control_idx_array`,
#'   `missing_idx_array`, `groups`, `units`, `times`.
#' @export
build_model_arrays <- function(df, groups, denominator_scale = 1e4,
                               allow_unbalanced_panel = FALSE) {
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
  filled <- array(FALSE, dim = c(K, D, N))

  k_idx <- match(df$group, groups)
  d_idx <- match(df$unit, units)
  n_idx <- match(df$time, times)
  cells <- cbind(k_idx, d_idx, n_idx)

  outcome_missing <- is.na(df$outcome)
  Y[cells[!outcome_missing, , drop = FALSE]] <-
    as.numeric(df$outcome[!outcome_missing])
  missing_idx[cells[outcome_missing, , drop = FALSE]] <- TRUE

  if ("denominator" %in% names(df)) {
    denom <- as.numeric(df$denominator)
    denom_ok <- !is.na(denom) & denom > 0
    denominators[cells[denom_ok, , drop = FALSE]] <-
      denom[denom_ok] / denominator_scale
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
    missing_idx <- missing_idx | !filled
    cli::cli_warn(
      "Unbalanced panel: {n_absent} of {K * D * N} cells are structurally
       absent (marked as missing)."
    )
  }

  list(
    Y = Y, denominators = denominators, control_idx_array = control_idx,
    missing_idx_array = missing_idx, groups = groups, units = units,
    times = times
  )
}

#' @export
print.bpnmf_data <- function(x, ...) {
  K <- length(x$groups)
  D <- length(x$units)
  N <- length(x$times)
  n_exposed <- sum(!x$control_idx_array)
  cli::cli_h1("bpnmf data ({x$type})")
  cli::cli_li("{K} group{?s}: {.val {x$groups}}")
  cli::cli_li("{D} unit{?s}, {N} time period{?s}")
  cli::cli_li("{n_exposed} exposed cell{?s}, {sum(x$missing_idx_array)} missing cell{?s}")
  invisible(x)
}
