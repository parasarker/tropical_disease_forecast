# 03_forecast_model.R
# iterative ensemble forecast for tropical disease (VL cases, São Paulo)
# corresponds to 4/6/2026 "Initial Iterative Ensemble Forecast" Milestone

#### outline ####
#   1. load target data
#   2. fit state-space model on all data up to today
#   3. get last latent state from the fit to use as the initial condition for next steps
#   4. forecast
#   5. format in EFI standard format
#   6. submit 

# NOTE on "starts from previous forecast":
#   This does NOT mean loading yesterday's output file.
#   It means each daily run re-fits the model with the latest data,
#   and uses X[T] (the last state of that fit) as the IC for the forecast.
#   Since our model is fast, we can re-fit from scratch each run.


#### step 0: load packages and data urls ####
library(rjags)
library(ecoforecastR)
library(tidyverse)
library(lubridate)
library(arrow)

model_id       <- "BU_visceral_leishmaniacs"   #name for leaderboard,, can change
forecast_date  <- Sys.Date()
n_ensemble     <- 100    # max ~100 for submission
# Forecast horizon (months): e.g. last obs Aug 2024 -> Sep 2024 ... Aug 2025
horizon        <- 12 # sept 2024 -> aug 2025

# data bucket URLS
s3_endpoint  <- "minio-s3.apps.shift.nerc.mghpcc.org"
target_url   <- "https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-read/challenges/project_id=bu4cast/targets/tropical-disease-targets.csv"
drivers_path <- "bu4cast-ci-read/challenges/project_id=bu4cast/drivers/"

# Monthly forecast met from stage2: TMP (monthly mean C) + APCP (mm/month total).
monthly_met_future <- function(raw, model_site_id) {
  raw <- dplyr::filter(raw, as.character(site_id) %in% as.character(model_site_id))
  tmp <- raw |>
    dplyr::filter(variable == "TMP") |>
    dplyr::mutate(month = as.Date(lubridate::floor_date(datetime, "month"))) |>
    dplyr::group_by(site_id, parameter, month) |>
    dplyr::summarise(TMP = mean(prediction), .groups = "drop") |>
    dplyr::group_by(parameter, month) |>
    dplyr::summarise(TMP = mean(TMP), .groups = "drop")
  apcp <- raw |>
    dplyr::filter(variable == "APCP") |>
    dplyr::mutate(month = as.Date(lubridate::floor_date(datetime, "month"))) |>
    dplyr::group_by(site_id, parameter, month) |>
    dplyr::summarise(APCP = sum(prediction), .groups = "drop") |>
    dplyr::group_by(parameter, month) |>
    dplyr::summarise(APCP = mean(APCP), .groups = "drop")
  out <- dplyr::full_join(tmp, apcp, by = c("parameter", "month"))
  dplyr::mutate(out, month = as.Date(as.character(month)))
}

# ------------------------------------------------------------------------------
#  THIS BLOCK SHOULD BE MADE INTO A SHARED FUNCTION BETWEEN 
# ------------------------------------------------------------------------------
#### step 1: load current target data ####
# get monthly disease, met and population data
source("00_clean_and_plot_monthly_inputs.R")
d <- monthly_data(start_date = as.Date("2007-01-01"), end_date = as.Date("2024-08-01"))
df <- d$monthly_merged

# pick one site from d$site_ids --> we should add the option to include all HERE
model_site_id <- d$site_ids[[1]]
model_site_monthly <- df %>%
  filter(site_id == model_site_id) %>%
  arrange(month)

## we should add the option to run sites one at a time or together
## maybe something like if(all_sites==TRUE){ RUN ALL SITES } else{ PICK }
##
##
##

y <- model_site_monthly$observation
time <- model_site_monthly$month
# Single temperature covariate comparable to monthly mean TMP on the forecast side
temp   <- (model_site_monthly$tmin_c + model_site_monthly$tmax_c) / 2
precip <- model_site_monthly$prec_mm
pop_scaled <- model_site_monthly$pop_scaled #scaled per 100k people
n_obs      <- length(y)

# First month after last observation through +horizon months (e.g. Sep 2024–Aug 2025 if data end Aug 2024)
forecast_start <- max(time) %m+% months(1)
future_months  <- seq(forecast_start, by = "month", length.out = horizon)
stopifnot(length(future_months) == horizon)

# ------------------------------------------------------------------------------

#### step 3: load future met drivers (stage2 from S3) ####
stage2 <- arrow::s3_bucket(paste0(drivers_path, "stage2/"), endpoint_override = s3_endpoint, anonymous = TRUE)
met_future <- arrow::open_dataset(stage2) |> dplyr::collect()
met_future <- met_future[as.character(met_future$site_id) == as.character(model_site_id), , drop = FALSE]

refs <- sort(unique(as.Date(as.character(met_future$reference_datetime))))
fd   <- forecast_start
reference_date <- if (length(refs) == 0L) {
  as.Date(NA)
} else if (any(refs <= fd)) {
  max(refs[refs <= fd])
} else {
  min(refs)
}
met_future <- met_future[as.Date(as.character(met_future$reference_datetime)) == reference_date, , drop = FALSE]

#### step 4: forecast-year population
pop_site <- d$pop_site |>
  dplyr::filter(site_id == model_site_id) |>
  dplyr::distinct(year, pop_est) |>
  dplyr::arrange(year)

last_year   <- max(pop_site$year)
last_pop    <- pop_site$pop_est[pop_site$year == last_year][1]
growth_rate <- mean(diff(pop_site$pop_est) / head(pop_site$pop_est, -1))

forecast_years <- (last_year + 1):(last_year + ceiling(horizon / 12) + 1)
pop_projected <- data.frame(
  year    = forecast_years,
  pop_est = last_pop * (1 + growth_rate)^(forecast_years - last_year)
)
pop_site <- rbind(pop_site, pop_projected)

#### step 5: fit model to current data ####

# data list for JAGS
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

# JAGS model string
StateSpace = "
model{
  for(t in 1:n_obs){
    y[t] ~ dnegbin(p[t], r)
    p[t] <- r / (r + mu[t])
    log(mu[t]) <- x[t]
  }
  for(t in 2:n_obs){
    x[t] ~ dnorm(mu_proc[t], tau_add)
    mu_proc[t] <- alpha + phi * x[t-1] +
                  beta_temp * temp[t] +
                  beta_precip * precip[t] +
                  beta_pop * pop[t]
  }
  x[1] ~ dnorm(x_ic, tau_ic)
  tau_add ~ dgamma(a_add, r_add)
  r ~ dgamma(1, 0.01)
  alpha ~ dnorm(0, 0.01)
  phi ~ dnorm(0, 0.01)
  beta_temp ~ dnorm(0, 0.01)
  beta_precip ~ dnorm(0, 0.01)
  beta_pop ~ dnorm(0, 0.01)}"

# Define initial model state
nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y, length(y), replace=TRUE)
  init[[i]] <- list(
    tau_add = 1 / max(stats::var(diff(log(y.samp + 1))), 1e-12),
    r = 10
  )
}

# Send info to JAGS
j.model <- jags.model(
  file = textConnection(StateSpace),
  data = jags_data,
  inits = init,
  n.chains = nchain)

# burn-in
jags.out <- coda.samples(
  model = j.model,
  variable.names = c("tau_add", "r"),
  n.iter = 2000)

# diagnostics
plot(jags.out)
# dic.samples(j.model, 2000)

# full sample 
jags.out <- coda.samples(
  model = j.model,
  variable.names = c(
    "x", "alpha", "phi",
    "beta_temp", "beta_precip", "beta_pop",
    "tau_add", "r"
  ),
  n.iter = 10000)

# extract posterior matrix
out <- as.matrix(jags.out)


#### step 6: get initial conditions from last fitted state ####

# grab all columns that correspond to latent state x
x.cols <- grep("^x\\[", colnames(out))

# notes!
# x[T] = last latent state = initial conditons for the forecast
# this is a distribution (one value per MCMC sample), not a single number !
x_T_all <- out[, max(x.cols)]

# subsample to n_ensemble members to propagate IC uncertainty into the forecast
idx    <- sample(nrow(out), n_ensemble, replace = nrow(out) < n_ensemble)
x_T    <- x_T_all[idx]       # IC vector, length n_ensemble
params <- out[idx, ]          # matching parameter samples for each ensemble member

#### step 7: run ensemble forecast forward ####
# future_months and horizon set in step 1

# get future met drivers (ensemble members x horizon)
met_monthly_forecast <- monthly_met_future(met_future, model_site_id)
n_met_members <- dplyr::n_distinct(met_monthly_forecast$parameter)
met_param <- (seq_len(n_ensemble) - 1) %% n_met_members
mo <- as.Date(format(future_months, "%Y-%m-01"))
mo_df <- tibble::tibble(month = mo)

temp_future   <- matrix(NA_real_, nrow = horizon, ncol = n_ensemble)
precip_future <- matrix(NA_real_, nrow = horizon, ncol = n_ensemble)

for (e in seq_len(n_ensemble)) {
  sub <- dplyr::filter(met_monthly_forecast, parameter == met_param[e]) |>
    dplyr::mutate(month = as.Date(as.character(month)))
  z <- dplyr::left_join(mo_df, sub, by = "month")
  temp_future[, e]   <- z$TMP
  precip_future[, e] <- z$APCP
}

# get future population
pop_future <- pop_site$pop_est[match(as.integer(format(future_months, "%Y")), pop_site$year)] / 1e6

# initialize forecast storage
x_fc <- matrix(NA, nrow = horizon, ncol = n_ensemble)  # latent state
y_fc <- matrix(NA, nrow = horizon, ncol = n_ensemble)  # predicted cases

# run forecast for each ensemble member
for (e in seq_len(n_ensemble)) {
  x_prev <- x_T[e]
  for (t in seq_len(horizon)) {
    mu_proc <- params[e, "alpha"] +
      params[e, "phi"] * x_prev +
      params[e, "beta_temp"] * temp_future[t, e] +
      params[e, "beta_precip"] * precip_future[t, e] +
      params[e, "beta_pop"] * pop_future[t]
    x_fc[t, e] <- stats::rnorm(1, mu_proc, sd = 1 / sqrt(params[e, "tau_add"]))
    x_prev <- x_fc[t, e]
  }
  
  # observation model -- convert latent state to predicted cases
  # now x_fc[t,e] is defined for every t
  for(t in 1:horizon){
    mu_t <- exp(x_fc[t, e])
    p_t  <- params[e, "r"] / (params[e, "r"] + mu_t)
    y_fc[t, e] <- rnbinom(1, size = params[e, "r"], prob = p_t)
  }
}


#### step 8: format in EFI standard format ####

# Required columns:
#   model_id, reference_datetime, datetime, site_id,
#   family ("ensemble"), parameter (1:n_ensemble),
#   variable ("cases"), prediction

# One row per [datetime x ensemble member]
# reference_datetime = reference_date (stage2 issue used; matches loaded drivers)
# datetime = sequence of future months

# forecast_long <- tibble(
#   model_id           = model_id,
#   reference_datetime = forecast_date,
#   datetime           = rep(future_months, each = n_ensemble),
#   site_id            = "SP",   # or aggregate site identifier
#   family             = "ensemble",
#   parameter          = rep(1:n_ensemble, times = horizon),
#   variable           = "cases",
#   prediction         = as.vector(t(y_fc))
# )

#### step 8b: visualize fitted values + forecast with CIs ####

# extract fitted latent states
x_fitted_cols <- matrix(out[idx, grep("^x\\[", colnames(out))], nrow = n_ensemble)
mu_fitted <- exp(x_fitted_cols)

# compute quantiles
fitted_ci  <- apply(mu_fitted, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

# build plot data frames
fitted_df <- tibble(
  month  = time,
  lo     = fitted_ci["2.5%",],
  median = fitted_ci["50%",],
  hi     = fitted_ci["97.5%",],
  obs    = y
)

# per lead time (works for horizon == 1; apply(..., 1, quantile) alone can return a vector)
forecast_df <- tibble(
  month  = future_months,
  lo     = apply(y_fc, 1, function(v) stats::quantile(v, 0.025, na.rm = TRUE)),
  median = apply(y_fc, 1, function(v) stats::quantile(v, 0.5, na.rm = TRUE)),
  hi     = apply(y_fc, 1, function(v) stats::quantile(v, 0.975, na.rm = TRUE))
)

# plot
# zoom to last year of historical + forecast
last_year_start <- max(time) %m-% months(12)

fitted_zoom <- filter(fitted_df, month >= last_year_start)

x_rng <- range(c(fitted_zoom$month, forecast_df$month), na.rm = TRUE)

ggplot() +
  geom_ribbon(data = fitted_zoom, aes(x = month, ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.3) +
  geom_line(data = fitted_zoom, aes(x = month, y = median), color = "steelblue", linewidth = 0.8) +
  geom_point(data = fitted_zoom, aes(x = month, y = obs), color = "black", size = 1.5) +
  geom_ribbon(data = forecast_df, aes(x = month, ymin = lo, ymax = hi), fill = "tomato", alpha = 0.3) +
  geom_line(data = forecast_df, aes(x = month, y = median), color = "tomato", linewidth = 1) +
  geom_vline(xintercept = as.numeric(max(time)), linetype = "dashed", color = "grey40") +
  scale_x_date(limits = x_rng, date_breaks = "1 month", date_labels = "%b %Y") +
  labs(title = "VL Forecast", x = "Month", y = "Cases") +
  theme_bw()

#### step 9: write forecast file and submit ####

# build forecast dataframe
forecast_long <- tibble(
  model_id           = model_id,
  reference_datetime = reference_date,
  datetime           = rep(future_months, times = n_ensemble),
  site_id            = model_site_id,
  family             = "ensemble",
  parameter          = rep(1:n_ensemble, each = horizon),
  variable           = "cases",
  prediction         = as.vector(y_fc)
)

head(forecast_long)

# filename convention: "<category>-<reference_datetime>-<model_id>.csv.gz"
forecast_file <- paste0("disease-", reference_date, "-", model_id, ".csv.gz")

# write compressed csv
write_csv(forecast_long, forecast_file)
file.exists(forecast_file) # did it work?

# submit to BU challenge S3 bucket 
# COMMENTED OUT UNTIL WE ARE READY!
# s3_submit <- arrow::s3_bucket("bu4cast-ci-write/challenges/project_id=bu4cast/forecasts/", endpoint_override = s3_endpoint, anonymous = FALSE)

# upload file
#aws.s3::put_object(file = forecast_file, object = forecast_file, bucket = "bu4cast-ci-write",region = "", base_url = s3_endpoint)
