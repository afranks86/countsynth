# Treatment-effect regression: lme4-style formula parsing and design
# construction for the flexible stage-2 / joint treatment block. A formula
# like  ~ 1 + event_time + (1 + event_time | group)  replaces the legacy
# group/unit/group:unit hierarchy with a linear fixed-effect surface plus
# random-effect terms whose per-level coefficients shrink toward that
# surface. Rows are the exposed cells in canonical flat order (utils.R), the
# same ordering that defines treatment_kt_z and cut provenance.

# Variables te_covariate_frame() always supplies, so a formula can use them
# without any covariate table. `time` is built too but is a Date and cannot
# enter a design; time_idx / event_time are its modelable forms.
TE_BUILTIN_VARS <- c("group", "unit", "time_idx", "event_time")

#' Treatment-effect regression options
#'
#' Configures the optional covariate model for the treatment effect used by
#' both joint inference and cut stage 2. When `formula` is `NULL` (the
#' default) the legacy hierarchy
#' `~ 0 + (1 | group) + (1 | unit) + (1 | group:unit)` is used via the
#' original model code path, so results are exactly backward compatible.
#'
#' @param formula RHS-only formula in lme4-style syntax, e.g.
#'   `~ 1 + event_time + (1 + event_time | group)`. Fixed-effect terms define
#'   the regression surface; each `(expr | g)` term adds per-level
#'   coefficients that shrink toward the surface. Available variables:
#'   `group`, `unit`, `time_idx` (1-based period index), `event_time`
#'   (periods since the unit's first treated period, 0 at onset), and any
#'   column of `covariates`. Grouping factors may be `group`, `unit`,
#'   interactions like `group:unit`, or any categorical covariate.
#'   An iid exposed-cell effect (the legacy `treatment_kt_z`, HalfNormal(0.1)
#'   scale) is always included and need not be written in the formula.
#' @param covariates Optional data frame (or CSV path) of cell-level
#'   covariates, joined onto the exposed cells by whichever of `unit`,
#'   `time`, `group` columns it contains (e.g. the post-ban change in
#'   distance to the nearest abortion clinic, keyed by `unit` or by
#'   `unit` + `time`). Covariates must be non-missing on every exposed cell
#'   they are used for.
#' @param standardize Center and scale continuous covariates before building
#'   the design (default `TRUE`); plots report functions on the original
#'   scale.
#' @param coef_prior_scale Prior sd for the fixed-effect coefficients,
#'   `te_beta ~ normal(0, coef_prior_scale)` (recycled to the number of
#'   design columns).
#' @param re_prior_scale Half-normal prior sd for each random-effect term's
#'   per-predictor scale (recycled across term predictors). Matches the
#'   legacy HalfNormal(1) scales at the default.
#' @return A `countsynth_te_opts` object.
#' @export
countsynth_te_opts <- function(formula = NULL, covariates = NULL,
                          standardize = TRUE, coef_prior_scale = 1,
                          re_prior_scale = 1) {
  if (!is.null(formula)) {
    if (is.character(formula)) {
      formula <- stats::as.formula(paste("~", sub("^\\s*~", "", formula)))
    }
    checkmate::assert_formula(formula)
    if (length(formula) != 2) {
      cli::cli_abort(
        "treatment_effects formula must be one-sided, e.g.
         {.code ~ 1 + event_time + (1 | group)}."
      )
    }
  }
  if (!is.null(covariates) && !is.character(covariates)) {
    checkmate::assert_data_frame(covariates, min.cols = 2)
  }
  checkmate::assert_flag(standardize)
  checkmate::assert_numeric(
    coef_prior_scale,
    lower = .Machine$double.xmin, min.len = 1, any.missing = FALSE
  )
  checkmate::assert_numeric(
    re_prior_scale,
    lower = .Machine$double.xmin, min.len = 1, any.missing = FALSE
  )
  if (is.null(formula) && !is.null(covariates)) {
    cli::cli_abort(
      "treatment_effects covariates were given without a {.field formula}."
    )
  }
  new_countsynth_class(
    list(
      formula = formula, covariates = covariates, standardize = standardize,
      coef_prior_scale = coef_prior_scale, re_prior_scale = re_prior_scale
    ),
    "countsynth_te_opts"
  )
}

# Extract every (expr | group) bar term from a formula RHS expression.
find_bars <- function(expr) {
  if (is.call(expr)) {
    op <- as.character(expr[[1]])
    if (op == "|") {
      return(list(expr))
    }
    if (op == "(") {
      return(find_bars(expr[[2]]))
    }
    if (op == "+" || op == "-") {
      bars <- find_bars(expr[[2]])
      if (length(expr) > 2) {
        bars <- c(bars, find_bars(expr[[3]]))
      }
      return(bars)
    }
  }
  list()
}

# Remove bar terms from a formula RHS expression; NULL when nothing is left.
strip_bars <- function(expr) {
  if (is.call(expr)) {
    op <- as.character(expr[[1]])
    if (op == "(" && is.call(expr[[2]]) &&
      identical(as.character(expr[[2]][[1]]), "|")) {
      return(NULL)
    }
    if ((op == "+" || op == "-") && length(expr) == 3) {
      lhs <- strip_bars(expr[[2]])
      rhs <- strip_bars(expr[[3]])
      if (is.null(rhs)) {
        return(lhs)
      }
      if (is.null(lhs)) {
        return(if (op == "-") call(op, rhs) else rhs)
      }
      return(call(op, lhs, rhs))
    }
  }
  expr
}

deparse1_ <- function(expr) paste(deparse(expr), collapse = " ")

#' Exposed-cell covariate frame
#'
#' One row per exposed cell in canonical flat order, holding the built-in
#' variables (`group`, `unit`, `time`, `time_idx`, `event_time`) plus any
#' user covariates joined by the key columns present in `covariates`.
#' @keywords internal
te_covariate_frame <- function(data, covariates = NULL) {
  D <- length(data$units)
  N <- length(data$times)
  control_flat <- flatten_kdn(data$control_idx_array)
  exp_cell <- which(!control_flat)
  sub <- kdn_from_flat(exp_cell, D, N)

  # First treated period per unit (any group); NA for never-treated units.
  exposed_any <- apply(!data$control_idx_array, c(2, 3), any)
  first_exposed <- apply(exposed_any, 1, function(r) {
    if (any(r)) which(r)[[1]] else NA_integer_
  })

  frame <- tibble::tibble(
    group = factor(data$groups[sub$k], levels = data$groups),
    unit = factor(data$units[sub$d], levels = data$units),
    time = data$times[sub$n],
    time_idx = sub$n,
    event_time = sub$n - first_exposed[sub$d]
  )

  if (!is.null(covariates)) {
    if (is.character(covariates)) {
      checkmate::assert_file_exists(covariates)
      covariates <- utils::read.csv(covariates, check.names = FALSE)
    }
    keys <- intersect(c("unit", "time", "group"), names(covariates))
    if (length(keys) == 0) {
      cli::cli_abort(
        "treatment_effects covariates must contain at least one key column
         among {.val {c('unit', 'time', 'group')}}."
      )
    }
    reserved <- setdiff(
      intersect(names(covariates), names(frame)), keys
    )
    if (length(reserved) > 0) {
      cli::cli_abort(
        "treatment_effects covariates use reserved column name{?s}
         {.val {reserved}}."
      )
    }
    if ("time" %in% keys) {
      covariates$time <- parse_time_column(covariates$time, "auto", "time")
    }
    if (anyDuplicated(covariates[keys])) {
      cli::cli_abort(
        "treatment_effects covariates have duplicate rows for key{?s}
         {.val {keys}}."
      )
    }
    join_frame <- frame
    join_frame$unit <- as.character(join_frame$unit)
    join_frame$group <- as.character(join_frame$group)
    joined <- dplyr::left_join(
      join_frame, covariates,
      by = keys
    )
    extra <- setdiff(names(covariates), keys)
    frame[extra] <- joined[extra]
  }
  frame
}

check_te_var_types <- function(frame, vars) {
  for (v in vars) {
    if (!v %in% names(frame)) {
      cli::cli_abort(
        "treatment_effects formula uses {.val {v}}, which is neither a
         built-in variable ({.val {TE_BUILTIN_VARS}}) nor a covariate
         column."
      )
    }
    col <- frame[[v]]
    if (inherits(col, "Date")) {
      cli::cli_abort(
        "treatment_effects formula uses the Date column {.val {v}};
         use {.field time_idx} or {.field event_time} (or a numeric
         covariate) instead."
      )
    }
    if (!is.numeric(col) && !is.factor(col) && !is.character(col) &&
      !is.logical(col)) {
      cli::cli_abort(
        "treatment_effects covariate {.val {v}} has unsupported type
         {.cls {class(col)}}."
      )
    }
    if (anyNA(col)) {
      cli::cli_abort(
        "treatment_effects covariate {.val {v}} is missing on
         {sum(is.na(col))} exposed cell{?s}."
      )
    }
  }
  invisible(NULL)
}

#' Build the treatment-effect regression design
#'
#' Parses the lme4-style formula in `opts`, evaluates it on the exposed-cell
#' covariate frame, and returns the fixed-effect design `X`, the per-term
#' random-effect designs, and the standardization metadata used to report
#' regression functions on the original covariate scale.
#'
#' @param opts A [countsynth_te_opts()] object with a non-`NULL` formula.
#' @param data A `countsynth_data` object.
#' @return A `countsynth_te_design` object.
#' @keywords internal
build_te_design <- function(opts, data) {
  checkmate::assert_class(opts, "countsynth_te_opts")
  if (is.null(opts$formula)) {
    return(NULL)
  }
  rhs <- opts$formula[[2]]
  bars <- find_bars(rhs)
  fixed_expr <- strip_bars(rhs)

  frame <- te_covariate_frame(data, opts$covariates)
  # Character covariates act as factors everywhere (design and grouping).
  for (nm in names(frame)) {
    if (is.character(frame[[nm]])) {
      frame[[nm]] <- factor(frame[[nm]])
    }
  }

  model_vars <- unique(c(
    all.vars(fixed_expr %||% quote(0)),
    unlist(lapply(bars, function(b) all.vars(b[[2]]))),
    unlist(lapply(bars, function(b) all.vars(b[[3]])))
  ))
  check_te_var_types(frame, model_vars)

  # Standardize continuous model variables (grouping-only factors are
  # untouched); centers/scales let plots map back to the original units.
  raw_frame <- frame
  centers <- numeric(0)
  scales <- numeric(0)
  if (opts$standardize) {
    for (v in model_vars) {
      col <- frame[[v]]
      # Indicators keep their 0/1 coding: a "per SD" rescaling of a dummy
      # reads worse than the raw contrast, and constants are left to the
      # rank check below.
      is_indicator <- is.numeric(col) && all(col %in% c(0, 1))
      if (is.numeric(col) && !is_indicator && stats::sd(col) > 0) {
        s <- stats::sd(col)
        m <- mean(col)
        frame[[v]] <- (col - m) / s
        centers[[v]] <- m
        scales[[v]] <- s
      }
    }
  }

  n <- nrow(frame)
  X <- if (is.null(fixed_expr)) {
    matrix(0, nrow = n, ncol = 0)
  } else {
    f <- stats::as.formula(call("~", fixed_expr), env = baseenv())
    stats::model.matrix(f, frame)
  }
  attr(X, "assign") <- NULL
  attr(X, "contrasts") <- NULL
  # A rank-deficient fixed design leaves directions the data cannot inform
  # (the classic case: a unit-level covariate that is constant because only
  # one unit is exposed). The prior would keep the sampler running while the
  # coefficients stayed unidentified, so refuse it here instead.
  if (ncol(X) > 0) {
    rank <- qr(X)$rank
    if (rank < ncol(X)) {
      constant <- setdiff(
        colnames(X)[apply(X, 2, function(cl) stats::sd(cl) == 0)],
        "(Intercept)"
      )
      cli::cli_abort(c(
        "treatment_effects fixed-effect design is rank deficient
         (rank {rank} < {ncol(X)} columns).",
        i = if (length(constant) > 0) {
          "Constant on the exposed cells: {.val {constant}}."
        } else {
          "Two or more design columns are collinear on the exposed cells."
        },
        i = "Only exposed cells enter this design, so unit-level covariates
             are constant unless several units are treated."
      ))
    }
  }

  terms <- lapply(bars, function(b) {
    g <- eval(b[[3]], envir = frame)
    if (is.character(g)) g <- factor(g)
    if (!is.factor(g)) {
      cli::cli_abort(
        "treatment_effects grouping factor {.code {deparse1_(b[[3]])}} must be
         categorical; got {.cls {class(g)}}."
      )
    }
    g <- droplevels(g)
    if (nlevels(g) == 1) {
      cli::cli_warn(c(
        "treatment_effects grouping factor {.code {deparse1_(b[[3]])}} has a
         single level on the exposed cells.",
        i = "Its coefficients are confounded with the corresponding fixed
             effects and only the prior separates them; drop the term or
             group by something that varies."
      ))
    }
    zf <- stats::as.formula(call("~", b[[2]]), env = baseenv())
    Z <- stats::model.matrix(zf, frame)
    attr(Z, "assign") <- NULL
    attr(Z, "contrasts") <- NULL
    if (ncol(Z) == 0) {
      cli::cli_abort(
        "treatment_effects term
         {.code ({deparse1_(b[[2]])} | {deparse1_(b[[3]])})} has no
         predictors."
      )
    }
    list(
      label = deparse1_(b[[3]]),
      levels = levels(g),
      index = as.integer(g),
      Z = Z,
      predictors = colnames(Z)
    )
  })

  new_countsynth_class(
    list(
      formula = opts$formula,
      # Parsed pieces are kept so the plotting layer can rebuild the same
      # design matrices on a prediction grid.
      fixed_expr = fixed_expr,
      bars = bars,
      standardize = opts$standardize,
      frame = frame,
      raw_frame = raw_frame,
      X = X,
      x_names = colnames(X),
      terms = terms,
      centers = centers,
      scales = scales,
      coef_prior_scale = rep_len(opts$coef_prior_scale, ncol(X)),
      re_prior_scale = rep_len(
        opts$re_prior_scale,
        sum(vapply(terms, function(t) ncol(t$Z), integer(1)))
      ),
      n_exposed = n
    ),
    "countsynth_te_design"
  )
}

#' @export
print.countsynth_te_design <- function(x, ...) {
  cli::cli_h1("countsynth treatment-effect design")
  cli::cli_li("formula: {.code {deparse1_(x$formula)}}")
  cli::cli_li("{x$n_exposed} exposed cell{?s}")
  cli::cli_li("fixed effects ({length(x$x_names)}): {.val {x$x_names}}")
  for (t in x$terms) {
    cli::cli_li(
      "({toString(t$predictors)} | {t$label}): {length(t$levels)} level{?s}"
    )
  }
  if (length(x$centers) > 0) {
    cli::cli_li("standardized: {.val {names(x$centers)}}")
  }
  invisible(x)
}

# Stan data fields for the treatment-effect regression block; the empty
# variant (te_reg = 0) selects the legacy hierarchy in the Stan models.
te_stan_fields <- function(design, n_exposed) {
  if (is.null(design)) {
    return(list(
      te_reg = 0L, P = 0L,
      X = matrix(0, nrow = n_exposed, ncol = 0),
      J = 0L, Qtot = 0L, Utot = 0L,
      L = integer(0), Q = integer(0),
      Z = matrix(0, nrow = n_exposed, ncol = 0),
      re_level = matrix(0L, nrow = 0, ncol = n_exposed),
      te_beta_prior_scale = numeric(0),
      te_re_prior_scale = numeric(0)
    ))
  }
  stopifnot(design$n_exposed == n_exposed)
  L <- vapply(design$terms, function(t) length(t$levels), integer(1))
  Q <- vapply(design$terms, function(t) ncol(t$Z), integer(1))
  Z <- do.call(cbind, c(
    lapply(design$terms, `[[`, "Z"),
    list(matrix(0, nrow = n_exposed, ncol = 0))
  ))
  re_level <- do.call(rbind, c(
    lapply(design$terms, function(t) matrix(t$index, nrow = 1)),
    list(matrix(0L, nrow = 0, ncol = n_exposed))
  ))
  list(
    te_reg = 1L, P = ncol(design$X),
    X = unname(design$X),
    J = length(design$terms), Qtot = sum(Q), Utot = sum(L * Q),
    L = as.array(L), Q = as.array(Q),
    Z = unname(Z),
    re_level = re_level,
    te_beta_prior_scale = as.array(design$coef_prior_scale),
    te_re_prior_scale = as.array(design$re_prior_scale)
  )
}
