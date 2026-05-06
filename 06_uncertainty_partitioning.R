# 06_uncertainty_partitioning.R
# Uncertainty partitioning for VL forecast
# This script builds the forecast 4 times, adding one uncertainty source at a time:
#     I   = Initial Conditions only
#     IP  = + Parameter uncertainty
#     IPD = + Driver uncertainty
#     IPDE= + Process error


# either run 02_historical_timeseries.R first and uncomment saveRDS at the bottom OR this script will re-fit the model for site 1 if no RDS is found

#### load packages ####
library(rjags)
library(ecoforecastR)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(lubridate)
library(cowplot)

source("00_clean_and_plot_monthly_inputs.R")

site_ids <- c(150210L, 500270L, 230440L, 170210L, 261110L, 355030L, 312770L, 310620L)
site_names <- c(
  "150210" = "Cametá (PA)",
  "500270" = "Campo Grande (MS)",
  "230440" = "Fortaleza (CE)",
  "170210" = "Araguaína (TO)",
  "261110" = "Petrolina (PE)",
  "355030" = "São Paulo (SP)",
  "312770" = "Governador Valadares (MG)",
  "310620" = "Belo Horizonte (MG)"
)

n_ensemble <- 1000
horizon    <- 12

dir.create("figures", showWarnings = FALSE)
dir.create("figures/uncertainty", showWarnings = FALSE)

#### JAGS model string ####
StateSpace <- "
model{
  for(t in 1:n_obs){
    y[t] ~ dnegbin(p[t], r)
    p[t] <- r / (r + mu[t])
    log(mu[t]) <- x[t]
  }
  for(t in 2:n_obs){
    x[t] ~ dnorm(mu_proc[t], tau_add)
    mu_proc[t] <- alpha + phi * x[t-1] +
                  beta_temp   * temp[t]   +
                  beta_precip * precip[t] +
                  beta_pop    * pop[t]
  }
  x[1]        ~ dnorm(x_ic, tau_ic)
  tau_add     ~ dgamma(a_add, r_add)
  r           ~ dgamma(1, 0.01)
  alpha       ~ dnorm(0, 0.01)
  phi         ~ dnorm(0, 0.01)
  beta_temp   ~ dnorm(0, 0.01)
  beta_precip ~ dnorm(0, 0.01)
  beta_pop    ~ dnorm(0, 0.01)
}"

#### load full dataset once ####
d  <- monthly_data(site_ids   = site_ids,
                   start_date = as.Date("2007-01-01"),
                   end_date   = as.Date("2024-08-01"))
df <- d$monthly_merged

#### helper: propagate ensemble forward ####
propagate <- function(x_ic_vec, alpha_vec, phi_vec, bt_vec, bp_vec,
                      bpop_vec, tau_vec, r_vec, temp_mat, precip_mat,
                      pop_last, add_process_error = FALSE) {
  ne   <- length(x_ic_vec)
  x_fc <- matrix(NA, nrow = horizon, ncol = ne)
  y_fc <- matrix(NA, nrow = horizon, ncol = ne)
  for (e in seq_len(ne)) {
    x_prev <- x_ic_vec[e]
    for (t in seq_len(horizon)) {
      mu_proc <- alpha_vec[e] + phi_vec[e] * x_prev +
        bt_vec[e]   * temp_mat[t, e] +
        bp_vec[e]   * precip_mat[t, e] +
        bpop_vec[e] * pop_last
      x_fc[t, e] <- if (add_process_error)
        rnorm(1, mu_proc, 1 / sqrt(tau_vec[e])) else mu_proc
      x_prev <- x_fc[t, e]
    }
    for (t in seq_len(horizon)) {
      mu_t       <- exp(x_fc[t, e])
      p_t        <- r_vec[e] / (r_vec[e] + mu_t)
      y_fc[t, e] <- rnbinom(1, size = r_vec[e], prob = p_t)
    }
  }
  list(x = x_fc, y = y_fc)
}

#### shared color palette ####
unc_colors <- c("IC" = "#2980B9", "Params" = "#E74C3C",
                "Drivers" = "#F39C12", "Process" = "#27AE60")
unc_labels <- c("IC" = "Initial Conditions", "Params" = "Parameters",
                "Drivers" = "Drivers (temp + precip + pop)", "Process" = "Process Error")

#### storage for combined panel ####
all_frac_plots <- list()

#### MAIN LOOP ####
for (site in site_ids) {
  
  site_label <- site_names[as.character(site)]
  message("\n=== Running uncertainty partitioning for: ", site_label, " ===")
  
  # ── prep site data ──────────────────────────────────────────────────────────
  model_site_monthly <- df %>% filter(site_id == site) %>% arrange(month)
  hasData <- !is.na(model_site_monthly$tmax_c) &
    !is.na(model_site_monthly$tmin_c)  &
    !is.na(model_site_monthly$prec_mm) &
    !is.na(model_site_monthly$pop_scaled)
  model_site_monthly <- model_site_monthly[hasData, ]
  
  if (nrow(model_site_monthly) < 12) {
    message("Skipping ", site_label, " — not enough data")
    next
  }
  
  y          <- model_site_monthly$observation
  temp       <- (model_site_monthly$tmin_c + model_site_monthly$tmax_c) / 2
  precip     <- model_site_monthly$prec_mm
  pop_scaled <- model_site_monthly$pop_scaled
  n_obs      <- length(y)
  pop_last   <- tail(pop_scaled[!is.na(pop_scaled)], 1)
  
  # ── fit or load MCMC ────────────────────────────────────────────────────────
  rds_file <- paste0("mcmc_", site, ".rds")
  
  if (file.exists(rds_file)) {
    message("Loading saved MCMC: ", rds_file)
    jags.out <- readRDS(rds_file)
  } else {
    message("Fitting JAGS model for ", site_label, "...")
    jags_data <- list(
      y = as.integer(y), n_obs = n_obs,
      temp = as.numeric(temp), precip = as.numeric(precip),
      pop = as.numeric(pop_scaled),
      x_ic = log(mean(y) + 1), tau_ic = 1, a_add = 1, r_add = 1
    )
    init <- lapply(1:3, function(i) {
      s <- sample(y, length(y), replace = TRUE)
      list(tau_add = 1 / max(var(diff(log(s + 1))), 1e-12), r = 10)
    })
    j.model <- jags.model(textConnection(StateSpace), data = jags_data,
                          inits = init, n.chains = 3, quiet = TRUE)
    coda.samples(j.model, c("tau_add", "r"), n.iter = 2000)
    jags.out <- coda.samples(j.model,
                             c("x", "alpha", "phi", "beta_temp",
                               "beta_precip", "beta_pop", "tau_add", "r"),
                             n.iter = 10000)
    saveRDS(jags.out, rds_file)
    message("Saved: ", rds_file)
  }
  
  out         <- as.matrix(jags.out)
  params_mean <- colMeans(out)
  
  # ── driver stats ─────────────────────────────────────────────────────────────
  temp_mean   <- mean(temp,   na.rm = TRUE)
  temp_sd     <- sd(temp,     na.rm = TRUE)
  precip_mean <- mean(precip, na.rm = TRUE)
  precip_sd   <- sd(precip,   na.rm = TRUE)
  
  # ── posterior draws ──────────────────────────────────────────────────────────
  idx         <- sample(nrow(out), n_ensemble, replace = TRUE)
  x_T_draws   <- out[idx, paste0("x[", n_obs, "]")]
  alpha_draws <- out[idx, "alpha"]
  phi_draws   <- out[idx, "phi"]
  bt_draws    <- out[idx, "beta_temp"]
  bp_draws    <- out[idx, "beta_precip"]
  bpop_draws  <- out[idx, "beta_pop"]
  tau_draws   <- out[idx, "tau_add"]
  r_draws     <- out[idx, "r"]
  
  temp_drv    <- matrix(rnorm(horizon * n_ensemble, temp_mean,   temp_sd),
                        nrow = horizon, ncol = n_ensemble)
  precip_drv  <- matrix(rnorm(horizon * n_ensemble, precip_mean, precip_sd),
                        nrow = horizon, ncol = n_ensemble)
  temp_clim   <- matrix(temp_mean,   nrow = horizon, ncol = n_ensemble)
  precip_clim <- matrix(precip_mean, nrow = horizon, ncol = n_ensemble)
  
  # ── 4 runs ───────────────────────────────────────────────────────────────────
  message("  Run 1: IC only...")
  N.I <- propagate(
    x_T_draws,
    rep(params_mean["alpha"], n_ensemble), rep(params_mean["phi"], n_ensemble),
    rep(params_mean["beta_temp"], n_ensemble), rep(params_mean["beta_precip"], n_ensemble),
    rep(params_mean["beta_pop"], n_ensemble), rep(params_mean["tau_add"], n_ensemble),
    rep(params_mean["r"], n_ensemble), temp_clim, precip_clim, pop_last, FALSE
  )$y
  
  message("  Run 2: IC + Params...")
  N.IP <- propagate(
    x_T_draws, alpha_draws, phi_draws, bt_draws, bp_draws,
    bpop_draws, tau_draws, r_draws, temp_clim, precip_clim, pop_last, FALSE
  )$y
  
  message("  Run 3: IC + Params + Drivers...")
  N.IPD <- propagate(
    x_T_draws, alpha_draws, phi_draws, bt_draws, bp_draws,
    bpop_draws, tau_draws, r_draws, temp_drv, precip_drv, pop_last, FALSE
  )$y
  
  message("  Run 4: IC + Params + Drivers + Process...")
  N.IPDE <- propagate(
    x_T_draws, alpha_draws, phi_draws, bt_draws, bp_draws,
    bpop_draws, tau_draws, r_draws, temp_drv, precip_drv, pop_last, TRUE
  )$y
  
  # ── variance partitioning ────────────────────────────────────────────────────
  var_I    <- apply(N.I,    1, var, na.rm = TRUE)
  var_IP   <- apply(N.IP,   1, var, na.rm = TRUE)
  var_IPD  <- apply(N.IPD,  1, var, na.rm = TRUE)
  var_IPDE <- apply(N.IPDE, 1, var, na.rm = TRUE)
  
  V_IC      <- var_I
  V_param   <- pmax(0, var_IP   - var_I)
  V_driver  <- pmax(0, var_IPD  - var_IP)
  V_process <- pmax(0, var_IPDE - var_IPD)
  V_total   <- var_IPDE
  
  frac_df <- data.frame(
    horizon = 1:horizon,
    IC      = V_IC      / V_total,
    Params  = V_param   / V_total,
    Drivers = V_driver  / V_total,
    Process = V_process / V_total
  ) %>%
    pivot_longer(-horizon, names_to = "source", values_to = "fraction") %>%
    mutate(source = factor(source, levels = c("IC", "Params", "Drivers", "Process")))
  
  abs_df <- data.frame(
    horizon = 1:horizon,
    IC      = V_IC,
    Params  = V_param,
    Drivers = V_driver,
    Process = V_process
  ) %>%
    pivot_longer(-horizon, names_to = "source", values_to = "variance") %>%
    mutate(source = factor(source, levels = c("IC", "Params", "Drivers", "Process")))
  
  # ── CI envelopes ─────────────────────────────────────────────────────────────
  ci_fn   <- function(mat) apply(mat, 1, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
  ci_I    <- ci_fn(N.I)
  ci_IP   <- ci_fn(N.IP)
  ci_IPD  <- ci_fn(N.IPD)
  ci_IPDE <- ci_fn(N.IPDE)
  h_seq   <- 1:horizon
  
  # ── per-site plots ───────────────────────────────────────────────────────────
  p_env <- ggplot() +
    geom_ribbon(aes(x = h_seq, ymin = ci_IPDE[1,], ymax = ci_IPDE[3,],
                    fill = "Process"), alpha = 0.4) +
    geom_ribbon(aes(x = h_seq, ymin = ci_IPD[1,],  ymax = ci_IPD[3,],
                    fill = "Drivers"), alpha = 0.5) +
    geom_ribbon(aes(x = h_seq, ymin = ci_IP[1,],   ymax = ci_IP[3,],
                    fill = "Params"),  alpha = 0.6) +
    geom_ribbon(aes(x = h_seq, ymin = ci_I[1,],    ymax = ci_I[3,],
                    fill = "IC"),      alpha = 0.8) +
    geom_line(aes(x = h_seq, y = ci_IPDE[2,]), color = "black", linewidth = 0.8) +
    scale_fill_manual(values = unc_colors, labels = unc_labels, name = NULL) +
    scale_x_continuous(breaks = 1:horizon) +
    labs(title = "Forecast Uncertainty Envelopes",
         x = "Forecast horizon (months)", y = "VL Cases (95% CI)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
  
  p_frac <- ggplot(frac_df, aes(x = horizon, y = fraction, fill = source)) +
    geom_area(alpha = 0.85, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = unc_colors, labels = unc_labels, name = NULL) +
    scale_x_continuous(breaks = 1:horizon) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "Relative Uncertainty Contributions",
         x = "Forecast horizon (months)", y = "Fraction of total variance") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
  
  p_abs <- ggplot(abs_df, aes(x = horizon, y = variance, color = source)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = unc_colors, labels = unc_labels, name = NULL) +
    scale_x_continuous(breaks = 1:horizon) +
    labs(title = "Absolute Variance by Source",
         x = "Forecast horizon (months)", y = "Variance") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9)) +
    guides(color = guide_legend(nrow = 2))
  
  site_panel <- p_env / (p_frac | p_abs) +
    plot_layout(heights = c(1.2, 1)) +
    plot_annotation(title = paste("Uncertainty Partitioning —", site_label))
  
  outfile <- paste0("figures/uncertainty/uncertainty_", site, ".png")
  ggsave(outfile, site_panel, width = 14, height = 10, dpi = 300)
  message("Saved: ", outfile)
  
  # store relative fraction plot for combined panel
  all_frac_plots[[as.character(site)]] <- p_frac +
    labs(title = site_label) +
    theme(plot.title = element_text(size = 9, face = "bold"))
  
  # cleanup memory after each site
  rm(jags.out, out, N.I, N.IP, N.IPD, N.IPDE)
  gc()
}

#### COMBINED PANEL — relative contributions all 8 sites ####
combined <- wrap_plots(all_frac_plots, ncol = 4) +
  plot_annotation(title = "Relative Uncertainty Contributions — All Sites")

# shared legend
legend_plot <- ggplot(
  data.frame(source = factor(names(unc_labels), levels = names(unc_labels)),
             x = 1, y = 1),
  aes(x, y, fill = source)
) +
  geom_col() +
  scale_fill_manual(values = unc_colors, labels = unc_labels, name = NULL) +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 11)) +
  guides(fill = guide_legend(nrow = 1))

legend_grob <- cowplot::get_legend(legend_plot)

final_combined <- cowplot::plot_grid(
  combined, legend_grob,
  ncol = 1, rel_heights = c(1, 0.05)
)

ggsave("figures/uncertainty_partitioning_all_sites.png",
       final_combined, width = 16, height = 9, dpi = 300)