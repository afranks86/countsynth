# Rank-shrinkage diagnostics. The frame builders are tested directly on
# synthetic draw arrays: the sampling-dependent part (pulling group_weight
# and eff_rank out of a CmdStan fit) is a two-line extractor, while the part
# that can be wrong -- sorting inside the draw rather than across draws -- is
# all here.

# [draw, group, component] array from a list of per-draw weight matrices.
profile_array <- function(draws) {
  K <- nrow(draws[[1]])
  R <- ncol(draws[[1]])
  out <- array(NA_real_, dim = c(length(draws), K, R))
  for (i in seq_along(draws)) {
    out[i, , ] <- draws[[i]]
  }
  out
}

test_that("the profile is sorted within each draw, not across draws", {
  # Two draws carrying the same weights in different component slots: the
  # component labels are not estimable (an NMF is invariant to permuting
  # them), so only the order statistics are, and both draws must contribute
  # the same values to the same position.
  w <- c(0.6, 0.3, 0.08, 0.02)
  gw <- profile_array(list(
    matrix(w, nrow = 1),
    matrix(w[c(3, 1, 4, 2)], nrow = 1),
    matrix(w[c(4, 3, 2, 1)], nrow = 1)
  ))
  f <- component_weight_frame(gw, groups = "total")
  expect_equal(f$component, 1:4)
  expect_equal(f$median, w)
  # Zero spread: every draw holds the same multiset, so a within-draw sort
  # leaves nothing for the bands to disagree about. Sorting the posterior
  # summaries instead would also pass here -- the next test separates them.
  expect_equal(f$lower_95, w)
  expect_equal(f$upper_95, w)
})

test_that("sorting per draw is not the same as sorting posterior means", {
  # Draw 1 loads component 1, draw 2 loads component 2, equally hard. The
  # per-component posterior means are both 0.5, so sorting summaries would
  # report a flat profile; sorting within the draw correctly reports one
  # dominant component and one unused.
  gw <- profile_array(list(
    matrix(c(0.9, 0.1), nrow = 1),
    matrix(c(0.1, 0.9), nrow = 1)
  ))
  f <- component_weight_frame(gw, groups = "total")
  expect_equal(f$median, c(0.9, 0.1))
})

test_that("per-group profiles stay separate and carry the group labels", {
  gw <- profile_array(list(
    matrix(c(0.5, 0.5, 0.9, 0.1), nrow = 2, byrow = TRUE),
    matrix(c(0.5, 0.5, 0.9, 0.1), nrow = 2, byrow = TRUE)
  ))
  f <- component_weight_frame(gw, groups = c("total", "medicaid"))
  expect_equal(f$group, rep(c("total", "medicaid"), each = 2))
  expect_equal(f$median, c(0.5, 0.5, 0.9, 0.1))
  # No group names available (a raw CmdStan fit): positional labels.
  expect_equal(
    unique(component_weight_frame(gw)$group), c("group 1", "group 2")
  )
})

test_that("eff_rank is summarized per group", {
  er <- matrix(c(2, 2.5, 1, 1.2), ncol = 2, byrow = TRUE)
  f <- eff_rank_frame(er, groups = c("a", "b"))
  expect_equal(f$group, c("a", "b"))
  expect_equal(f$median, c(1.5, 1.85))
})

test_that("a fit without rank shrinkage is reported, not plotted", {
  no_shrink <- structure(
    list(stan_data = list(rank_shrink = 0L)), class = c("bpnmf_fit", "list")
  )
  expect_null(rank_shrinkage_source(no_shrink))
  expect_error(
    bpnmf_component_weight_summary(no_shrink), "not run with rank shrinkage"
  )
  # A cut fit is unwrapped to stage 1, which is where the factor block lives.
  cut_fit <- structure(
    list(stage1 = no_shrink), class = c("bpnmf_cut_fit", "list")
  )
  expect_null(rank_shrinkage_source(cut_fit))
})
