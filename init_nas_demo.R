# Purpose: Size-structured demographic scenario analysis — evaluates
#   matrix vs iter_bins equilibrium initialization across movement,
#   mortality, and growth scenarios. Same to the age-structured version.


# setup -------------------------------------------------------------------

library(here)
library(tidyverse)
library(cowplot)
library(patchwork)


# helper functions --------------------------------------------------------

# Growth transition matrix for a single region (gamma increments, terminal bin accumulates)
build_Xr <- function(n_sizes, mean_growth, cv_growth = 0.5) {
  Xr <- matrix(0, n_sizes, n_sizes)
  for (i in 1:(n_sizes - 1)) {
    mean_inc <- mean_growth * (1 - 0.7 * (i - 1) / (n_sizes - 1))
    sd_inc   <- mean_inc * cv_growth
    shape    <- (mean_inc / sd_inc)^2
    rate_g   <- mean_inc / sd_inc^2
    probs    <- rep(0, n_sizes)
    for (j in i:n_sizes) {
      probs[j] <- pgamma(j - i + 1, shape = shape, rate = rate_g) -
        pgamma(j - i,     shape = shape, rate = rate_g)
    }
    probs[n_sizes] <- probs[n_sizes] + (1 - sum(probs))  # remainder into terminal bin
    Xr[, i] <- probs
  }
  Xr[n_sizes, n_sizes] <- 1.0  # terminal bin stays
  return(Xr)
}

# Block-diagonal growth matrix across regions: X_full is (n_sizes*n_regions) x (n_sizes*n_regions)
build_X_full <- function(n_regions, n_sizes, X_list) {
  X_full <- matrix(0, n_sizes * n_regions, n_sizes * n_regions)
  for (r in 1:n_regions) {
    idx <- ((r - 1) * n_sizes + 1):(r * n_sizes)
    X_full[idx, idx] <- X_list[[r]]
  }
  X_full
}

# Full movement matrix T_full with psi[to, from] structure for each size bin.
# State vector is stacked [R1 sizes | R2 sizes | ...].
build_T_full <- function(n_regions, n_sizes, psi_list) {
  T_full <- matrix(0, n_sizes * n_regions, n_sizes * n_regions)
  for (r in 1:n_regions) {
    for (rp in 1:n_regions) {
      row_idx   <- ((r  - 1) * n_sizes + 1):(r  * n_sizes)
      col_idx   <- ((rp - 1) * n_sizes + 1):(rp * n_sizes)
      diag_vals <- sapply(1:n_sizes, function(l) psi_list[[l]][r, rp])
      T_full[row_idx, col_idx] <- diag(diag_vals)
    }
  }
  T_full
}

# Full transition matrix O = X_full %*% S_full %*% T_full
# natmort_mat: n_sizes x n_regions (instantaneous mortality rates)
build_O <- function(n_regions, n_sizes, X_full, natmort_mat, T_full) {
  S_surv <- exp(-natmort_mat)         # n_sizes x n_regions survival probabilities
  S_full <- diag(as.vector(S_surv))   # diagonal, stacked [R1; R2; ...]
  X_full %*% S_full %*% T_full
}

# Build recruitment entry vector r_tilde (recruits enter first size bin per region)
build_r_tilde <- function(n_sizes, n_regions, r0r) {
  r_tilde <- rep(0, n_sizes * n_regions)
  for (r in 1:n_regions) r_tilde[(r - 1) * n_sizes + 1] <- r0r[r]
  r_tilde
}

# matrix equilibrium: n* = (I - O)^{-1} r_tilde
get_equil_matrix <- function(O, r_tilde, n_sizes, n_regions) {
  I_mat <- diag(n_sizes * n_regions)
  n_eq  <- solve(I_mat - O, r_tilde)
  matrix(n_eq, n_sizes, n_regions)
}

# iter_bins equilibrium: n_{t+1} = O n_t + r_tilde until convergence
get_equil_iter_bins <- function(O, r_tilde, n_sizes, n_regions, n_iter = 500) {
  n_v <- rep(0, n_sizes * n_regions)
  for (t in 1:n_iter) n_v <- O %*% n_v + r_tilde
  matrix(n_v, n_sizes, n_regions)
}

# No-movement matrix: uses identity movement, otherwise same O structure.
# Analogous to no_move_all_* methods in the age-structured case.
get_equil_no_move_all <- function(n_regions, n_sizes, X_full, natmort_mat, r_tilde) {
  T_id  <- diag(n_sizes * n_regions)
  O_nm  <- build_O(n_regions, n_sizes, X_full, natmort_mat, T_id)
  get_equil_matrix(O_nm, r_tilde, n_sizes, n_regions)
}

# setup scenarios ---------------------------------------------------------

n_regions <- 2
n_sizes   <- 15

# Growth scenarios (region-specific mean growth rate, same CV)
growth_list <- list(build_Xr(n_sizes, 2, 0.5), build_Xr(n_sizes, 2, 0.5))  # both fast

# Natural mortality scenarios (n_sizes x n_regions matrices of rates)
natmort_list <- list(
  eq   = cbind(rep(0.25, n_sizes), rep(0.25, n_sizes)),
  uneq = cbind(rep(0.20, n_sizes), rep(0.25, n_sizes))
)

# Movement scenarios
# psi_list: list of n_sizes matrices, each n_regions x n_regions with psi[to, from]
move_list <- list(

  # No movement
  no = lapply(1:n_sizes, function(l) diag(1, n_regions)),

  # Random/equal mixing (same rate for all size bins)
  rand = lapply(1:n_sizes, function(l) matrix(1 / n_regions, n_regions, n_regions)),

  # Symmetric exchange (30% move between regions, size-invariant)
  oneway = {
    rate <- 0.3
    lapply(1:n_sizes, function(l)
      matrix(c(1 - rate, rate,
               rate,     1 - rate), nrow = n_regions, byrow = TRUE))
  },

  # Ontogenetic movement: logistic increase in emigration from R1 to R2 with size
  onto = {
    midpoint    <- n_sizes / 2
    steepness   <- 0.1
    max_rate    <- 0.8
    leave_rates <- max_rate / (1 + exp(-steepness * (1:n_sizes - midpoint)))
    lapply(1:n_sizes, function(l)
      matrix(c(1 - leave_rates[l], leave_rates[l],
               1 - max_rate,       max_rate), nrow = n_regions, byrow = TRUE))
  }
)

# Recruitment distribution scenarios
r0r_list <- list(
  eq  = c(0.5, 0.5),
  src = c(0.7, 0.3)   # source-sink (R1 dominant recruiter)
)

# Methods to compare (analogous to init_age_strc in age-structured version)
init_methods <- c("matrix", "iter_bins", "no_move_all")


# factorial design --------------------------------------------------------

factorial_scenarios <- expand.grid(
  movement = names(move_list),
  natmort  = names(natmort_list),
  r0r      = names(r0r_list),
  stringsAsFactors = FALSE
)

scenario_list <- vector("list", nrow(factorial_scenarios))
for (i in 1:nrow(factorial_scenarios)) {
  nm <- paste(factorial_scenarios$movement[i],
              factorial_scenarios$natmort[i],
              factorial_scenarios$r0r[i], sep = "_")
  scenario_list[[i]] <- list(
    name     = nm,
    movement = move_list[[factorial_scenarios$movement[i]]],
    natmort  = natmort_list[[factorial_scenarios$natmort[i]]],
    r0r      = r0r_list[[factorial_scenarios$r0r[i]]]
  )
}
names(scenario_list) <- sapply(scenario_list, function(x) x$name)


  # run scenarios -----------------------------------------------------------

init_nss_results <- data.frame()

init_methods <- c("matrix", "iter_bins", "no_move_all")

for (i in seq_along(scenario_list)) {

  tmp_scenario <- scenario_list[[i]]

  X_full  <- build_X_full(n_regions, n_sizes, growth_list)
  T_full  <- build_T_full(n_regions, n_sizes, tmp_scenario$movement)
  O       <- build_O(n_regions, n_sizes, X_full, tmp_scenario$natmort, T_full)
  r_tilde <- build_r_tilde(n_sizes, n_regions, tmp_scenario$r0r)

  spec_rad <- max(abs(eigen(O, only.values = TRUE)$values))
  if (spec_rad >= 1) warning(paste("Spectral radius >= 1 for scenario:", tmp_scenario$name))

  for (method in init_methods) {
    n_eq <- switch(method,
                   matrix   = get_equil_matrix(O, r_tilde, n_sizes, n_regions),
                   iter_bins    = get_equil_iter_bins(O, r_tilde, n_sizes, n_regions),
                   no_move_all      = get_equil_no_move_all(n_regions, n_sizes, X_full, tmp_scenario$natmort, r_tilde)
    )

    tmp_df <- reshape2::melt(n_eq) %>%
      rename(Size = Var1, Region = Var2) %>%
      mutate(scenario = tmp_scenario$name, method = method)
    init_nss_results <- rbind(tmp_df, init_nss_results)
  }
}


# tidy results ------------------------------------------------------------

init_nss_results <- init_nss_results %>%
  separate(scenario, into = c("movement", "natmort", "r0r"),
           sep = "_", remove = FALSE) %>%
  mutate(
    method   = factor(method, levels = c("matrix", "iter_bins", "no_move_all")),
    movement = factor(movement, levels = c("no", "rand", "oneway", "onto"),
                      labels = paste("move:", c("no", "rand", "oneway", "onto"))),
    natmort  = factor(natmort, levels = c("eq", "uneq"),
                      labels = paste("natmort:", c("eq", "uneq")))
  )

# Relative difference vs matrix (reference method)
init_nss_results <- init_nss_results %>%
  group_by(scenario, Region, Size) %>%
  mutate(
    ref_val  = value[method == "matrix"],
    rel_diff = (value - ref_val) / ref_val
  ) %>%
  ungroup()


# plot: region 1 ----------------------------------------------------------

# coloring and line types
cb_pal <- ggthemes::colorblind_pal()(8)  # full 8-color palette
nas_colors <- cb_pal[c(1, 3, 4)]  # skip position 2
nas_ltys   <- c("solid", "dashed", "dotted")  # match NAA lty order for these 3


main_r1 <- init_nss_results %>%
  filter(Region == 1, r0r == "src") %>%
  ggplot(aes(x = Size, y = value, color = method, lty = method)) +
  geom_line(lwd = 1.3) +
  ggh4x::facet_grid2(movement ~ natmort, scales = "free_y", independent = "y") +
  labs(color = "Method", lty = "Method", y = "Equilibrium Abundance (Region 1)") +
  scale_color_manual(values = setNames(nas_colors, c("no_move_all", "iter_bins", "matrix"))) +
  scale_fill_manual(values  = setNames(nas_colors, c("no_move_all", "iter_bins", "matrix"))) +
  scale_linetype_manual(values = setNames(nas_ltys, c("no_move_all", "iter_bins", "matrix"))) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")

make_inset <- function(region, move_label, natmort_label) {
  init_nss_results %>%
    filter(Region == region, r0r == "src", Size == n_sizes,
           movement == move_label, natmort == natmort_label) %>%
    ggplot(aes(x = Size, y = value, fill = method)) +
    geom_col(position = position_dodge(), color = "black", alpha = 0.85, lwd = 0.1) +
    scale_color_manual(values = setNames(nas_colors, c("no_move_all", "iter_bins", "matrix"))) +
    scale_fill_manual(values  = setNames(nas_colors, c("no_move_all", "iter_bins", "matrix"))) +
    scale_linetype_manual(values = setNames(nas_ltys, c("no_move_all", "iter_bins", "matrix"))) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

region_1 <- main_r1 +
  inset_element(make_inset(1, "move: no",     "natmort: eq"),   left = 0.28, bottom = 0.85, right = 0.43, top = 0.99) +
  inset_element(make_inset(1, "move: no",     "natmort: uneq"), left = 0.81, bottom = 0.85, right = 0.96, top = 0.99) +
  inset_element(make_inset(1, "move: rand",   "natmort: eq"),   left = 0.28, bottom = 0.59, right = 0.43, top = 0.73) +
  inset_element(make_inset(1, "move: rand",   "natmort: uneq"), left = 0.81, bottom = 0.59, right = 0.96, top = 0.73) +
  inset_element(make_inset(1, "move: oneway", "natmort: eq"),   left = 0.28, bottom = 0.33, right = 0.43, top = 0.47) +
  inset_element(make_inset(1, "move: oneway", "natmort: uneq"), left = 0.81, bottom = 0.33, right = 0.96, top = 0.47) +
  inset_element(make_inset(1, "move: onto",   "natmort: eq"),   left = 0.28, bottom = 0.08, right = 0.43, top = 0.22) +
  inset_element(make_inset(1, "move: onto",   "natmort: uneq"), left = 0.81, bottom = 0.08, right = 0.96, top = 0.22)


# plot: region 2 ----------------------------------------------------------

main_r2 <- init_nss_results %>%
  filter(Region == 2, r0r == "src") %>%
  ggplot(aes(x = Size, y = value, color = method, lty = method)) +
  geom_line(lwd = 1.3) +
  ggh4x::facet_grid2(movement ~ natmort, scales = "free_y", independent = "y") +
  labs(color = "Method", lty = "Method", y = "Equilibrium Abundance (Region 2)") +
  scale_color_manual(values = setNames(nas_colors, c("no_move_all", "iter_bins", "matrix"))) +
  scale_fill_manual(values  = setNames(nas_colors, c("no_move_all", "iter_bins", "matrix"))) +
  scale_linetype_manual(values = setNames(nas_ltys, c("no_move_all", "iter_bins", "matrix"))) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")

region_2 <- main_r2 +
  inset_element(make_inset(2, "move: no",     "natmort: eq"),   left = 0.28, bottom = 0.85, right = 0.43, top = 0.99) +
  inset_element(make_inset(2, "move: no",     "natmort: uneq"), left = 0.81, bottom = 0.85, right = 0.96, top = 0.99) +
  inset_element(make_inset(2, "move: rand",   "natmort: eq"),   left = 0.28, bottom = 0.59, right = 0.43, top = 0.73) +
  inset_element(make_inset(2, "move: rand",   "natmort: uneq"), left = 0.81, bottom = 0.59, right = 0.96, top = 0.73) +
  inset_element(make_inset(2, "move: oneway", "natmort: eq"),   left = 0.28, bottom = 0.33, right = 0.43, top = 0.47) +
  inset_element(make_inset(2, "move: oneway", "natmort: uneq"), left = 0.81, bottom = 0.33, right = 0.96, top = 0.47) +
  inset_element(make_inset(2, "move: onto",   "natmort: eq"),   left = 0.28, bottom = 0.08, right = 0.43, top = 0.22) +
  inset_element(make_inset(2, "move: onto",   "natmort: uneq"), left = 0.81, bottom = 0.08, right = 0.96, top = 0.22)

legend <- cowplot::get_legend(main_r1 + theme(legend.position = "top") + labs(fill = "Method"))
comb   <- cowplot::plot_grid(region_1, region_2, labels = c("A", "B"), label_size = 30)
comb2  <- cowplot::plot_grid(legend, comb, ncol = 1, rel_heights = c(0.05, 0.95))

ggsave(
  here("size_struct_demo_plot.png"),
  comb2,
  width = 16, height = 12
)
