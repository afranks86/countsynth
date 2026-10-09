# Rank-shrinkage diagnostics: the two figures that answer "was this rank
# large enough?" without a rank sweep. Both read the factor block, which in
# cut mode lives in stage 1 -- stage 2 never sees it -- so a cut fit is
# unwrapped before anything is read.
#
# Every quantity plotted here is a permutation-invariant function of
# group_weight. Component identity is not estimable: an NMF is invariant to
# relabeling its components and nothing in the prior orders them, so
# group_weight[k, 3] means a different curve in each chain (the same reason
# DEFAULT_GATE_PARAMS leaves the factor parameters ungated). What is
# estimable is the profile's order statistics -- sorted inside each draw,
# never across draws -- and its inverse Simpson index.

# The fit holding the shrinkage parameters, or NULL when the fit did not use
# them.
rank_shrinkage_source <- function(x) {
  if (inherits(x, "countsynth_cut_fit")) {
    x <- x$stage1
  }
  if (inherits(x, "countsynth_fit")) {
    if (!isTRUE(x$stan_data$rank_shrink == 1L)) {
      return(NULL)
    }
    return(list(
      fit = x$fit, groups = x$data$groups, stan_data = x$stan_data
    ))
  }
  # A raw CmdStan fit: the parameter's presence is the only signal available.
  vars <- tryCatch(
    x$metadata()$stan_variables,
    error = function(e) character()
  )
  if (!"group_weight" %in% vars) {
    return(NULL)
  }
  list(fit = x, groups = NULL, stan_data = NULL)
}

rank_source_or_abort <- function(x) {
  rank_shrinkage_source(x) %||% cli::cli_abort(c(
    "This fit was not run with rank shrinkage.",
    i = "Set {.field model.rank_shrinkage} (see
         {.fn countsynth_rank_shrinkage_opts}) and refit."
  ))
}

# Group labels, falling back to positional names for a raw CmdStan fit.
rank_group_labels <- function(groups, K) {
  if (length(groups) == K) {
    as.character(groups)
  } else {
    sprintf("group %d", seq_len(K))
  }
}

# Median and nested 67% / 95% bands, one row per column of `m`.
rank_posterior_band <- function(m) {
  q <- function(p) apply(m, 2, stats::quantile, probs = p, names = FALSE)
  tibble::tibble(
    median = apply(m, 2, stats::median),
    lower_67 = q(0.165), upper_67 = q(0.835),
    lower_95 = q(0.025), upper_95 = q(0.975)
  )
}

# Profile draws as a [draw, group, component] array and the effective rank as
# [draw, group]; read together so both describe the same draws.
rank_draws <- function(src) {
  rv <- posterior::as_draws_rvars(
    src$fit$draws(variables = c("group_weight", "eff_rank"))
  )
  list(
    profile = posterior::draws_of(rv$group_weight),
    eff_rank = posterior::draws_of(rv$eff_rank)
  )
}

# Order statistics of the shared profile. Sorting within the draw is what
# makes this a posterior over something: "the largest weight" is well
# defined, "component 3's weight" is not.
component_weight_frame <- function(profile, groups = NULL) {
  d <- dim(profile)
  if (length(d) != 3L) {
    cli::cli_abort("Expected a [draw, group, component] array of draws.")
  }
  n_draws <- d[[1]]
  labels <- rank_group_labels(groups, d[[2]])
  dplyr::bind_rows(lapply(seq_len(d[[2]]), function(k) {
    m <- matrix(profile[, k, ], nrow = n_draws, ncol = d[[3]])
    # apply() returns [component, draw]; byrow refills it as [draw, rank].
    sorted <- matrix(
      apply(m, 1, sort, decreasing = TRUE),
      nrow = n_draws, ncol = d[[3]], byrow = TRUE
    )
    dplyr::bind_cols(
      tibble::tibble(group = labels[[k]], component = seq_len(d[[3]])),
      rank_posterior_band(sorted)
    )
  }))
}

eff_rank_frame <- function(eff_rank, groups = NULL) {
  m <- as.matrix(eff_rank)
  dplyr::bind_cols(
    tibble::tibble(group = rank_group_labels(groups, ncol(m))),
    rank_posterior_band(m)
  )
}

#' Posterior of the shared component-popularity profile
#'
#' The order statistics of `group_weight[k]` -- the shared profile of the
#' finite-HDP weight prior ([countsynth_rank_shrinkage_opts()]) -- sorted within
#' each posterior draw. One row per (group, position), with the position's
#' posterior median and 67% / 95% credible bounds.
#'
#' @param x A `countsynth_fit` / `countsynth_cut_fit` fit with `rank_shrinkage` set.
#' @return A tibble with `group`, `component`, `median`, `lower_67`,
#'   `upper_67`, `lower_95`, `upper_95`.
#' @export
countsynth_component_weight_summary <- function(x) {
  src <- rank_source_or_abort(x)
  component_weight_frame(rank_draws(src)$profile, src$groups)
}

#' Posterior of the effective number of factor components in use
#'
#' `eff_rank[k] = 1 / sum_r group_weight[k, r]^2`, the inverse Simpson index
#' of the shared profile: it reads `m` when `m` components split the mass
#' evenly and the remaining `R - m` are unused. A posterior well below the
#' fitted rank says the truncation was generous enough; one pressing against
#' it says to raise the rank and refit.
#'
#' @inheritParams countsynth_component_weight_summary
#' @return A tibble with `group`, `median`, `lower_67`, `upper_67`,
#'   `lower_95`, `upper_95`.
#' @export
countsynth_eff_rank_summary <- function(x) {
  src <- rank_source_or_abort(x)
  eff_rank_frame(rank_draws(src)$eff_rank, src$groups)
}

#' Scree plot of the shared component-popularity profile
#'
#' The posterior of the sorted profile weights, one panel per group, with the
#' effective rank marked. Read it as a Bayesian scree plot: positions whose
#' interval has collapsed toward zero are components the group is not using,
#' and their presence is what makes the fit insensitive to the rank. If the
#' last position still holds appreciable weight, the rank was too small.
#'
#' The dotted line is `1 / R`, the weight every position would carry under a
#' flat profile that uses all `R` components.
#'
#' @inheritParams countsynth_component_weight_summary
#' @return A ggplot object.
#' @export
countsynth_component_weight_plot <- function(x) {
  src <- rank_source_or_abort(x)
  dr <- rank_draws(src)
  prof <- component_weight_frame(dr$profile, src$groups)
  eff <- eff_rank_frame(dr$eff_rank, src$groups)
  R <- max(prof$component)

  # The effective rank belongs in the strip label rather than a legend: it is
  # one number per panel and reads with the curve it summarizes.
  strip <- stats::setNames(
    sprintf(
      "%s \u2014 eff. rank %.1f [%.1f, %.1f]",
      eff$group, eff$median, eff$lower_95, eff$upper_95
    ),
    eff$group
  )
  prof$facet <- factor(strip[prof$group], levels = unname(strip))
  eff$facet <- factor(strip[eff$group], levels = unname(strip))

  ggplot2::ggplot(prof, ggplot2::aes(x = .data$component)) +
    ggplot2::geom_hline(
      yintercept = 1 / R, linetype = "dotted", color = "grey40"
    ) +
    ggplot2::geom_vline(
      data = eff, ggplot2::aes(xintercept = .data$median),
      linetype = "dashed", color = "#C44E52"
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = .data$lower_95, ymax = .data$upper_95),
      linewidth = 0.7, alpha = 0.4
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = .data$lower_67, ymax = .data$upper_67),
      linewidth = 1.8, alpha = 0.9
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$median),
      shape = 21, fill = "white", size = 2.4, stroke = 0.9
    ) +
    ggplot2::scale_x_continuous(
      breaks = if (R <= 15) seq_len(R) else ggplot2::waiver()
    ) +
    ggplot2::facet_wrap("facet") +
    ggplot2::labs(
      title = "Shared component-popularity profile",
      subtitle = paste(
        "Profile weights sorted within each draw; thick 67% CI, thin 95% CI.",
        sprintf(
          "Dotted: 1/R = %.3f (flat profile). Dashed: effective rank.",
          1 / R
        )
      ),
      x = "Position in the sorted profile", y = "Profile weight"
    ) +
    theme_countsynth()
}

# Facet labels for the two hyperparameters: the parameter name plus what
# moving it does, so the figure is readable without the help page.
RANK_HYPER_LABELS <- c(
  unit_weight_sd =
    "unit_weight_sd \u2014 how far a group's units depart from its profile",
  group_weight_mass =
    "group_weight_mass \u2014 DP mass per group (sparser when small)"
)

# Prior density functions as passed to Stan, or NULL for a fit whose data
# list is not available (a raw CmdStan fit). The two are different families:
# a Gamma on the DP mass, a half-normal on the log-scale spread.
rank_hyper_priors <- function(stan_data) {
  if (is.null(stan_data)) {
    return(NULL)
  }
  list(
    unit_weight_sd = function(x) {
      2 * stats::dnorm(x, 0, stan_data$unit_sd_scale)
    },
    group_weight_mass = function(x) {
      stats::dgamma(x, stan_data$group_mass_shape, stan_data$group_mass_rate)
    }
  )
}

#' Posterior of the two rank-shrinkage hyperparameters
#'
#' Densities for `unit_weight_sd` (how far each unit's component loadings
#' depart from its group's profile, multiplicatively) and
#' `group_weight_mass` (the DP mass behind the stick-breaking profile, which
#' sets how sparse it is -- one per group, drawn separately), with their
#' priors overlaid.
#'
#' The prior overlay is the point of the figure: these are what the shrinkage
#' estimates rather than assumes, so a posterior that has not moved off its
#' prior says the panel had little to say about how much structure the units
#' share -- and that the shrinkage is being driven by the hyperprior rather
#' than the data.
#'
#' Read the mass panel across groups as well as against the prior. Each group
#' gets its own density, and one sitting well away from the rest marks a group
#' whose profile is unlike the others' -- often one whose cells are mostly
#' censored, carrying little information about how many components it needs.
#'
#' @inheritParams countsynth_component_weight_summary
#' @param prior Overlay the prior densities (default `TRUE`).
#' @return A ggplot object.
#' @export
countsynth_rank_hyper_plot <- function(x, prior = TRUE) {
  src <- rank_source_or_abort(x)
  rv <- posterior::as_draws_rvars(
    src$fit$draws(variables = c("unit_weight_sd", "group_weight_mass"))
  )
  # Both hyperparameters are one value per group, so their draws carry a group
  # dimension. Keeping the groups apart is the point: one group's density
  # sitting away from the others is what says its profile, or how far its
  # units sit from it, is unlike the rest -- and pooling the draws into a
  # single density would hide exactly that.
  df <- dplyr::bind_rows(lapply(names(RANK_HYPER_LABELS), function(nm) {
    d <- posterior::draws_of(rv[[nm]])
    n_col <- if (length(dim(d)) > 1L) dim(d)[[2]] else 1L
    labels <- if (n_col == length(src$groups)) {
      src$groups
    } else if (n_col > 1L) {
      as.character(seq_len(n_col))
    } else {
      NA_character_
    }
    tibble::tibble(
      parameter = nm,
      group = rep(labels, each = dim(d)[[1]]),
      value = as.vector(d)
    )
  }))
  df$facet <- factor(
    RANK_HYPER_LABELS[df$parameter], levels = unname(RANK_HYPER_LABELS)
  )

  priors <- if (prior) rank_hyper_priors(src$stan_data) else NULL
  prior_df <- NULL
  if (!is.null(priors)) {
    prior_df <- dplyr::bind_rows(lapply(names(priors), function(nm) {
      v <- df$value[df$parameter == nm]
      # Over the posterior's own range: the prior is here for comparison, and
      # letting its tail set the axis would squeeze the posterior flat.
      grid <- seq(0, max(v), length.out = 200)
      tibble::tibble(
        parameter = nm, value = grid, density = priors[[nm]](grid)
      )
    }))
    prior_df$facet <- factor(
      RANK_HYPER_LABELS[prior_df$parameter],
      levels = unname(RANK_HYPER_LABELS)
    )
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_density(
      ggplot2::aes(fill = .data$group, color = .data$group),
      alpha = 0.35
    ) +
    ggplot2::scale_fill_discrete(na.value = "#4C72B0", name = NULL) +
    ggplot2::scale_color_discrete(na.value = "#4C72B0", name = NULL)
  if (!is.null(prior_df)) {
    p <- p +
      ggplot2::geom_line(
        data = prior_df, ggplot2::aes(y = .data$density),
        linetype = "dashed", color = "grey35", linewidth = 0.7
      )
  }
  p +
    ggplot2::facet_wrap("facet", scales = "free") +
    ggplot2::labs(
      title = "Rank-shrinkage hyperparameters",
      subtitle = if (is.null(prior_df)) {
        "Posterior densities."
      } else {
        "Solid: posterior. Dashed: prior."
      },
      x = "Value", y = "Density"
    ) +
    theme_countsynth()
}

# Standard figure set for the report. Names become the PNG filenames; they
# are written into figs/ppc/ because they answer the same question the
# PPC-driven rank sweep was answering.
rank_report_figures <- function(x) {
  list(
    rank_component_weight = countsynth_component_weight_plot(x),
    rank_hyperparameters = countsynth_rank_hyper_plot(x)
  )
}
