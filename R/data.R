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

validate_denominators <- function(df, outcomes, unit_col) {
  for (o in outcomes) {
    denom_col <- o$denominator_col
    if (!is.null(denom_col) && denom_col %in% names(df)) {
      vals <- df[[denom_col]]
      bad <- is.na(vals) | vals <= 0
      if (any(bad)) {
        bad_units <- unique(df[[unit_col]][bad])
        cli::cli_abort(
          "NaN or non-positive denominator in column {.field {denom_col}}:
           {sum(bad)} row{?s} affected (units: {.val {bad_units}})."
        )
      }
    }
  }
  invisible(NULL)
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

aggregate_temporal <- function(df, period) {
  months_per_period <- switch(period,
    monthly = 1L, bimonthly = 2L, quarterly = 3L, yearly = 12L,
    cli::cli_abort(
      "Unknown aggregation period {.val {period}}; use monthly, bimonthly,
       quarterly, or yearly."
    )
  )
  month <- as.integer(format(df$time, "%m"))
  period_code <- if (period == "yearly") {
    rep(1L, nrow(df))
  } else {
    (month - 1L) %/% months_per_period + 1L
  }
  df$.year <- as.integer(format(df$time, "%Y"))
  df$.period <- period_code
  has_denom <- "denominator" %in% names(df)

  df_agg <- df |>
    dplyr::group_by(.data$unit, .data$.year, .data$.period, .data$group) |>
    dplyr::summarise(
      outcome = sum(.data$outcome),
      treatment = max(.data$treatment),
      denominator = if (has_denom) mean(.data$denominator) else 1,
      .groups = "drop"
    )

  first_month <- (df_agg$.period - 1L) * months_per_period + 1L
  df_agg$time <- as.Date(
    sprintf("%d-%02d-01", df_agg$.year, first_month)
  )
  df_agg$start_date <- df_agg$time
  df_agg$end_date <- add_months(df_agg$time, months_per_period) - 1L

  df_agg$.year <- NULL
  df_agg$.period <- NULL
  dplyr::arrange(df_agg, .data$unit, .data$time, .data$group)
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
  validate_denominators(df, outcomes, schema$unit_col)
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

  if (config$aggregation$enabled) {
    df_long <- aggregate_temporal(df_long, config$aggregation$period)
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
