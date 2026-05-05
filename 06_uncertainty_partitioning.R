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

source("00_clean_and_plot_monthly_inputs.R")

site_ids <- c(150210L, 500270L, 230440L, 170210L, 261110L, 355030L, 312770L, 310620L)
site_names <- c(
  "150210" = "Camet\u00e1 (PA)",
  "500270" = "Campo Grande (MS)",
  "230440" = "Fortaleza (CE)",
  "170210" = "Aragua\u00edna (TO)",
  "261110" = "Petrolina (PE)",
  "355030" = "S\u00e3o Paulo (SP)",
  "312770" = "Governador Valadares (MG)",
  "310620" = "Belo Horizonte (MG)"
)

n_ensemble <- 1000   # more members = smoother variance estimates
horizon    <- 12     # months ahead — longer horizon shows partitioning better

dir.create("figures", showWarnings = FALSE)


#### load or refit MCMC ####
if (file.exists("mcmc_output.rds")) {
  message("Loading saved MCMC output...")
  jags.out <- readRDS("mcmc_output.rds")
  
  # also need the data to get met drivers and last observed state
  d  <- monthly_data(site_ids   = site_ids,
                     start_date = as.Date("2007-01-01"),
                     end_date   = as.Date("2024-08-01"))
  df <- d$monthly_merged
  
  model_site_id      <- site_ids[1]
  model_site_monthly <- df %>% filter(site_id == model_site_id) %>% arrange(month)
  hasData <- !is.na(model_site_monthly$tmax_c) &
    !is.na(model_site_monthly$tmin_c)  &
    !is.na(model_site_monthly$prec_mm) &
    !is.na(model_site_monthly$pop_scaled)
  model_site_monthly <- model_site_monthly[hasData, ]
  y          <- model_site_monthly$observation
  time       <- model_site_monthly$month
  temp       <- (model_site_monthly$tmin_c + model_site_monthly$tmax_c) / 2
  precip     <- model_site_monthly$prec_mm
  pop_scaled <- model_site_monthly$pop_scaled
  
} else {
  message("No mcmc_output.rds found — re-fitting model for site 1...")
  
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
  
  d  <- monthly_data(site_ids   = site_ids,
                     start_date = as.Date("2007-01-01"),
                     end_date   = as.Date("2024-08-01"))
  df <- d$monthly_merged
  
  model_site_id      <- site_ids[1]
  model_site_monthly <- df %>% filter(site_id == model_site_id) %>% arrange(month)
  hasData <- !is.na(model_site_monthly$tmax_c) &
    !is.na(model_site_monthly$tmin_c)  &
    !is.na(model_site_monthly$prec_mm) &
    !is.na(model_site_monthly$pop_scaled)
  model_site_monthly <- model_site_monthly[hasData, ]
  y          <- model_site_monthly$observation
  time       <- model_site_monthly$month
  temp       <- (model_site_monthly$tmin_c + model_site_monthly$tmax_c) / 2
  precip     <- model_site_monthly$prec_mm
  pop_scaled <- model_site_monthly$pop_scaled
  n_obs      <- length(y)
  
  jags_data <- list(
    y       = as.integer(y),
    n_obs   = n_obs,
    temp    = as.numeric(temp),
    precip  = as.numeric(precip),
    pop     = as.numeric(pop_scaled),
    x_ic    = log(mean(y) + 1),
    tau_ic  = 1,
    a_add   = 1,
    r_add   = 1
  )
  
  init <- lapply(1:3, function(i) {
    s <- sample(y, length(y), replace = TRUE)
    list(tau_add = 1 / max(var(diff(log(s + 1))), 1e-12), r = 10)
  })
  
  j.model <- jags.model(textConnection(StateSpace), data = jags_data,
                        inits = init, n.chains = 3, quiet = TRUE)
  coda.samples(j.model, c("tau_add", "r"), n.iter = 2000)
  jags.out <- coda.samples(j.model,
                           c("x", "alpha", "phi",
                             "beta_temp", "beta_precip", "beta_pop",
                             "tau_add", "r"),
                           n.iter = 10000)
  saveRDS(jags.out, "mcmc_output.rds")
  message("MCMC saved to mcmc_output.rds")
}

out    <- as.matrix(jags.out)
x.cols <- grep("^x\\[", colnames(out))
n_obs  <- length(y)


#### posterior means (used for deterministic runs) ####
params_mean <- colMeans(out)

# posterior mean of last latent state (IC)
x_T_mean <- params_mean[paste0("x[", n_obs, "]")]

# posterior SD of last latent state (IC uncertainty)
x_T_sd <- sd(out[, paste0("x[", n_obs, "]")])

# met driver climatology mean & sd for driver uncertainty
temp_mean   <- mean(temp,   na.rm = TRUE)
temp_sd     <- sd(temp,     na.rm = TRUE)
precip_mean <- mean(precip, na.rm = TRUE)
precip_sd   <- sd(precip,   na.rm = TRUE)
pop_last    <- tail(pop_scaled[!is.na(pop_scaled)], 1)


#### helper: propagate one ensemble forward h steps ####
propagate <- function(
    x_ic_vec,       # vector length n_ensemble: initial condition draws
    alpha_vec,      # vector: alpha draws (or rep of mean)
    phi_vec,        # vector: phi draws
    bt_vec,         # vector: beta_temp draws
    bp_vec,         # vector: beta_precip draws
    bpop_vec,       # vector: beta_pop draws
    tau_vec,        # vector: tau_add draws (process error precision)
    r_vec,          # vector: r (NegBin dispersion) draws
    temp_mat,       # matrix horizon x n_ensemble: driver draws
    precip_mat,     # matrix horizon x n_ensemble
    add_process_error = FALSE
) {
  ne     <- length(x_ic_vec)
  x_fc   <- matrix(NA, nrow = horizon, ncol = ne)
  y_fc   <- matrix(NA, nrow = horizon, ncol = ne)
  
  for (e in seq_len(ne)) {
    x_prev <- x_ic_vec[e]
    for (t in seq_len(horizon)) {
      mu_proc <- alpha_vec[e] +
        phi_vec[e]   * x_prev        +
        bt_vec[e]    * temp_mat[t, e]   +
        bp_vec[e]    * precip_mat[t, e] +
        bpop_vec[e]  * pop_last
      
      if (add_process_error) {
        x_fc[t, e] <- rnorm(1, mu_proc, 1 / sqrt(tau_vec[e]))
      } else {
        x_fc[t, e] <- mu_proc
      }
      x_prev <- x_fc[t, e]
    }
    # observation model: NegBin
    for (t in seq_len(horizon)) {
      mu_t       <- exp(x_fc[t, e])
      p_t        <- r_vec[e] / (r_vec[e] + mu_t)
      y_fc[t, e] <- rnbinom(1, size = r_vec[e], prob = p_t)
    }
  }
  list(x = x_fc, y = y_fc)
}


#### sample draws from posterior ####
idx <- sample(nrow(out), n_ensemble, replace = TRUE)

# IC: last latent state
x_T_draws <- out[idx, paste0("x[", n_obs, "]")]

# parameters
alpha_draws    <- out[idx, "alpha"]
phi_draws      <- out[idx, "phi"]
bt_draws       <- out[idx, "beta_temp"]
bp_draws       <- out[idx, "beta_precip"]
bpop_draws     <- out[idx, "beta_pop"]
tau_draws      <- out[idx, "tau_add"]
r_draws        <- out[idx, "r"]

# driver uncertainty: resample from historical climate distribution
temp_drv   <- matrix(rnorm(horizon * n_ensemble, temp_mean,   temp_sd),
                     nrow = horizon, ncol = n_ensemble)
precip_drv <- matrix(rnorm(horizon * n_ensemble, precip_mean, precip_sd),
                     nrow = horizon, ncol = n_ensemble)

# climatological drivers (mean only, no uncertainty)
temp_clim   <- matrix(temp_mean,   nrow = horizon, ncol = n_ensemble)
precip_clim <- matrix(precip_mean, nrow = horizon, ncol = n_ensemble)


#### RUN 1: Initial Conditions only ####
# fix params at posterior mean, fix drivers at mean, no process error
message("Run 1: IC only...")
N.I <- propagate(
  x_ic_vec  = x_T_draws,               # IC uncertainty ON
  alpha_vec  = rep(params_mean["alpha"],        n_ensemble),
  phi_vec    = rep(params_mean["phi"],          n_ensemble),
  bt_vec     = rep(params_mean["beta_temp"],    n_ensemble),
  bp_vec     = rep(params_mean["beta_precip"],  n_ensemble),
  bpop_vec   = rep(params_mean["beta_pop"],     n_ensemble),
  tau_vec    = rep(params_mean["tau_add"],      n_ensemble),
  r_vec      = rep(params_mean["r"],            n_ensemble),
  temp_mat   = temp_clim,
  precip_mat = precip_clim,
  add_process_error = FALSE
)$y

#### RUN 2: IC + Parameters ####
message("Run 2: IC + Parameters...")
N.IP <- propagate(
  x_ic_vec  = x_T_draws,               # IC uncertainty ON
  alpha_vec  = alpha_draws,            # param uncertainty ON
  phi_vec    = phi_draws,
  bt_vec     = bt_draws,
  bp_vec     = bp_draws,
  bpop_vec   = bpop_draws,
  tau_vec    = tau_draws,
  r_vec      = r_draws,
  temp_mat   = temp_clim,
  precip_mat = precip_clim,
  add_process_error = FALSE
)$y

#### RUN 3: IC + Parameters + Drivers ####
message("Run 3: IC + Parameters + Drivers...")
N.IPD <- propagate(
  x_ic_vec  = x_T_draws,
  alpha_vec  = alpha_draws,
  phi_vec    = phi_draws,
  bt_vec     = bt_draws,
  bp_vec     = bp_draws,
  bpop_vec   = bpop_draws,
  tau_vec    = tau_draws,
  r_vec      = r_draws,
  temp_mat   = temp_drv,              # driver uncertainty ON
  precip_mat = precip_drv,
  add_process_error = FALSE
)$y

#### RUN 4: IC + Parameters + Drivers + Process Error ####
message("Run 4: IC + Parameters + Drivers + Process Error...")
N.IPDE <- propagate(
  x_ic_vec  = x_T_draws,
  alpha_vec  = alpha_draws,
  phi_vec    = phi_draws,
  bt_vec     = bt_draws,
  bp_vec     = bp_draws,
  bpop_vec   = bpop_draws,
  tau_vec    = tau_draws,
  r_vec      = r_draws,
  temp_mat   = temp_drv,
  precip_mat = precip_drv,
  add_process_error = TRUE            # process error ON
)$y


#### compute variances at each horizon step ####
var_I    <- apply(N.I,    1, var, na.rm = TRUE)
var_IP   <- apply(N.IP,   1, var, na.rm = TRUE)
var_IPD  <- apply(N.IPD,  1, var, na.rm = TRUE)
var_IPDE <- apply(N.IPDE, 1, var, na.rm = TRUE)

# incremental contributions (each source = difference from previous run)
V_IC      <- var_I
V_param   <- pmax(0, var_IP   - var_I)
V_driver  <- pmax(0, var_IPD  - var_IP)
V_process <- pmax(0, var_IPDE - var_IPD)
V_total   <- var_IPDE

# relative contributions (fraction of total variance)
frac_df <- data.frame(
  horizon  = 1:horizon,
  IC       = V_IC      / V_total,
  Params   = V_param   / V_total,
  Drivers  = V_driver  / V_total,
  Process  = V_process / V_total
) %>%
  pivot_longer(-horizon, names_to = "source", values_to = "fraction") %>%
  mutate(source = factor(source,
                         levels = c("IC", "Params", "Drivers", "Process")))


# shared color palette
unc_colors <- c(
  "IC"      = "#2980B9",
  "Params"  = "#E74C3C",
  "Drivers" = "#F39C12",
  "Process" = "#27AE60"
)
unc_labels <- c(
  "IC"      = "Initial Conditions",
  "Params"  = "Parameters",
  "Drivers" = "Drivers (temp + precip)",
  "Process" = "Process Error"
)

#### PLOT 1: stacked CI envelopes (Activity 11 style) ####
# Draw from outermost (full) to innermost (IC only) so inner bands are visible
ci_fn <- function(mat) apply(mat, 1, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)

ci_I    <- ci_fn(N.I)
ci_IP   <- ci_fn(N.IP)
ci_IPD  <- ci_fn(N.IPD)
ci_IPDE <- ci_fn(N.IPDE)

h_seq <- 1:horizon
med   <- ci_IPDE[2, ]   # median from full run for reference line

p_env <- ggplot() +
  # outermost first so inner bands paint on top
  geom_ribbon(aes(x = h_seq, ymin = ci_IPDE[1,], ymax = ci_IPDE[3,],
                  fill = "Process"), alpha = 0.4) +
  geom_ribbon(aes(x = h_seq, ymin = ci_IPD[1,],  ymax = ci_IPD[3,],
                  fill = "Drivers"), alpha = 0.5) +
  geom_ribbon(aes(x = h_seq, ymin = ci_IP[1,],   ymax = ci_IP[3,],
                  fill = "Params"),  alpha = 0.6) +
  geom_ribbon(aes(x = h_seq, ymin = ci_I[1,],    ymax = ci_I[3,],
                  fill = "IC"),      alpha = 0.8) +
  geom_line(aes(x = h_seq, y = med), color = "black", linewidth = 0.8) +
  scale_fill_manual(values = unc_colors, labels = unc_labels,
                    name = "Uncertainty source") +
  scale_x_continuous(breaks = 1:horizon) +
  labs(
    title = "Forecast Uncertainty Envelopes",
    x = "Forecast horizon (months)",
    y = "VL Cases (95% CI)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")


#### PLOT 2: stacked area — relative variance contributions ####
p_frac <- ggplot(frac_df, aes(x = horizon, y = fraction, fill = source)) +
  geom_area(alpha = 0.85, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = unc_colors, labels = unc_labels,
                    name = "Uncertainty source") +
  scale_x_continuous(breaks = 1:horizon) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Relative Uncertainty Contributions",
    x = "Forecast horizon (months)",
    y = "Fraction of total variance"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.text  = element_text(size = 10)) +
  guides(fill = guide_legend(nrow = 2))


#### PLOT 3: absolute variance by source ####
abs_df <- data.frame(
  horizon = 1:horizon,
  IC      = V_IC,
  Params  = V_param,
  Drivers = V_driver,
  Process = V_process
) %>%
  pivot_longer(-horizon, names_to = "source", values_to = "variance") %>%
  mutate(source = factor(source, levels = c("IC", "Params", "Drivers", "Process")))

p_abs <- ggplot(abs_df, aes(x = horizon, y = variance, color = source)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = unc_colors, labels = unc_labels,
                     name = "Source") +
  scale_x_continuous(breaks = 1:horizon) +
  labs(
    title = "Absolute Variance by Uncertainty Source",
    x = "Forecast horizon (months)",
    y = "Variance"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")


#### combine and save ####
full_panel <- p_env / (p_frac | p_abs) +
  plot_layout(heights = c(1.2, 1))

ggsave("figures/uncertainty_partitioning.png", full_panel,
       width = 14, height = 12, dpi = 300)

print(full_panel)

