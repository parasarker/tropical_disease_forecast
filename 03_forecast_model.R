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
library(daymetr)

model_id       <- "BU_visceral_leishmaniacs"   #name for leaderboard,, can change
forecast_date  <- Sys.Date()
n_ensemble     <- 100    # max ~100 for submission 
horizon        <- 2      # months forward to forecast (adjust depending)

# data bucket URLS
s3_endpoint  <- "minio-s3.apps.shift.nerc.mghpcc.org"
target_url   <- "https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-read/challenges/project_id=bu4cast/targets/tropical-disease-targets.csv"
drivers_path <- "bu4cast-ci-read/challenges/project_id=bu4cast/drivers/"

# parquet met (dirs from repo data/stage3_temp, data/stage2_temp)
stage3_dir <- file.path("data", "stage3_temp")
stage2_dir <- file.path("data", "stage2_temp")

# convert met from kelvin to C
kelvin_to_celsius <- function(x) { as.numeric(x)-273.15 }

# monthly covariates from stage3 (one row per parameter x month)
monthly_stage3 <- function(s3, sp_sites) {
  s3 <- dplyr::filter(s3, .data$site_id %in% sp_sites)
  air <- s3 |>
    dplyr::filter(.data$variable == "air_temperature") |>
    # aggregate to monthly
    dplyr::mutate(month = as.Date(lubridate::floor_date(.data$datetime, "month"))) |>
    dplyr::group_by(.data$site_id, .data$parameter, .data$month) |>
    dplyr::summarise(
      # we want tMax and tMin because that's what stage2 has
      tMax = max(.data$prediction, na.rm = TRUE),
      tMin = min(.data$prediction, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      tMax = kelvin_to_celsius(.data$tMax), 
      tMin = kelvin_to_celsius(.data$tMin)
    )
  prec <- s3 |>
    dplyr::filter(.data$variable == "precipitation_flux") |>
    dplyr::mutate(month = as.Date(lubridate::floor_date(.data$datetime, "month"))) |>
    dplyr::group_by(.data$site_id, .data$parameter, .data$month) |>
    dplyr::summarise(
      # stage3 has precipitation flux, stage2 has accumulated precipitation (apcp)
      # these are both amount/month so we can probably treat as same
      precip = sum(.data$prediction, na.rm = TRUE), .groups = "drop")
  dplyr::full_join(air, prec, by = c("site_id", "parameter", "month")) |>
    dplyr::group_by(.data$parameter, .data$month) |>
    dplyr::summarise(
      tMax = mean(.data$tMax, na.rm = TRUE),
      tMin = mean(.data$tMin, na.rm = TRUE),
      precip = mean(.data$precip, na.rm = TRUE),
      .groups = "drop"
    )
}

# monthly tMax, tMin, precipitation from stage2  (one row per parameter x month).
monthly_stage2 <- function(s2, sp_sites) {
  s2 <- dplyr::filter(s2, .data$site_id %in% sp_sites)
  tmp <- s2 |>
    dplyr::filter(.data$variable %in% c("TMAX", "TMIN")) |>
    dplyr::mutate(month = as.Date(lubridate::floor_date(.data$datetime, "month"))) |>
    dplyr::group_by(.data$site_id, .data$parameter, .data$month, .data$variable) |>
    dplyr::summarise(prediction = mean(.data$prediction, na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(.data$parameter, .data$month, .data$variable) |>
    dplyr::summarise(prediction = mean(.data$prediction, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = .data$variable, values_from = .data$prediction)
  prec <- s2 |>
    dplyr::filter(.data$variable == "APCP") |>
    dplyr::mutate(month = as.Date(lubridate::floor_date(.data$datetime, "month"))) |>
    dplyr::group_by(.data$site_id, .data$parameter, .data$month) |>
    dplyr::summarise(APCP = sum(.data$prediction, na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(.data$parameter, .data$month) |>
    dplyr::summarise(APCP = mean(.data$APCP, na.rm = TRUE), .groups = "drop")
  out <- dplyr::full_join(tmp, prec, by = c("parameter", "month"))
  for (cn in c("TMAX", "TMIN", "APCP")) {
    if (!cn %in% names(out)) out[[cn]] <- NA_real_
  }
  # stage2 tmax/min are already in C
  out
}

#### step 1: load current target data ####

# pull latest disease observations
disease_targets <- read.csv(target_url)

# subset to Sao Paulo (site_id 350000–360000)
sp <- subset(disease_targets, site_id >= 350000 & site_id < 360000)
sp$month <- as.Date(format(as.Date(sp$datetime), "%Y-%m-01"))
sp_monthly <- aggregate(observation ~ month, data = sp, sum, na.rm = TRUE)

time <- sp_monthly$month
y    <- sp_monthly$observation

# get site_ids to filter met data later on
sp_sites <- unique(sp$site_id)


#### step 2: load historical met drivers  ####

daymet0 <- download_daymet(site = "SaoPaulo", start = 2024, end = 2025, internal = TRUE)
daymet <- daymet0$data

# make date col
daymet$date <- as.Date(paste(daymet$year, daymet$yday, sep = "-"), "%Y-%j")
daymet$month <- as.Date(format(daymet$date, "%Y-%m-01"))

#aggregate to monthly to match disease data
tMin_monthly <- aggregate(daymet$tmin..deg.c., by = list(daymet$month), mean, na.rm = TRUE)
tMax_monthly <- aggregate(daymet$tmax..deg.c., by = list(daymet$month), mean, na.rm = TRUE)
precip_monthly <- aggregate(daymet$prcp..mm.day., by = list(daymet$month), sum, na.rm = TRUE)

time_month <- as.Date(format(time, "%Y-%m-01"))

# extract temp and precip data
tMin <- tMin_monthly[,2][match(time_month, tMin_monthly[,1])]
tMax <- tMax_monthly[,2][match(time_month, tMax_monthly[,1])]
precip <- precip_monthly[,2][match(time_month, precip_monthly[,1])]

# was going to use stage3 data (monthly, 30 ensembles) but it only has data from the past two months right now, so I have used the same daymet data that we used in the previous script as a placeholder right now...
# heres the code I wouldve used to bring that data in tho 
# s3 <- arrow::s3_bucket(paste0(drivers_path, "stage3/"), endpoint_override = s3_endpoint, anonymous = TRUE)
# met_hist <- arrow::open_dataset(s3) |> filter(site_id %in% sp_sites) |> collect()

# time_month <- as.Date(format(time, "%Y-%m-01"))
# 
# ds3 <- arrow::open_dataset(stage3_dir) |>
#   dplyr::filter(site_id %in% sp_sites) |>
#   dplyr::collect()
# 
# met_s3 <- monthly_stage3(ds3, sp_sites)
# met_fit <- dplyr::filter(met_s3, .data$parameter == 0)
# tMax   <- met_fit$tMax[match(time_month, met_fit$month)]
# tMin   <- met_fit$tMin[match(time_month, met_fit$month)]
# precip <- met_fit$precip[match(time_month, met_fit$month)]
# 
# # remove time points where climate data is unavailable
# hasData <- !is.na(tMax) & !is.na(tMin) & !is.na(precip)
# y      <- y[hasData]
# tMax   <- tMax[hasData]
# tMin   <- tMin[hasData]
# precip <- precip[hasData]
# time   <- time[hasData]
# n      <- length(y)

# met_hist <- arrow::open_dataset(s3) |> filter(site_id %in% sp_sites) |> collect()

# Align met drivers to disease time vector, handle NAs


#### step 3: load future met drivers (stage2) ####

# stage2 <- arrow::s3_bucket(paste0(drivers_path, "stage2/"),endpoint_override = s3_endpoint,anonymous = TRUE)
# met_future <- arrow::open_dataset(s3_stage2) |> filter(site_id %in% sp_sites) |> collect()

#### step 4: load population drivers ####
pop_df <- read.csv("data/state_pop_by_year.csv")
pop_sp <- subset(pop_df, state == "SP")
pop_sp <- pop_sp[order(pop_sp$year), ]

# extend population forward for forecast horizon using recent growth rate
last_year   <- tail(pop_sp$year, 1)
last_pop    <- tail(pop_sp$pop_est, 1)
growth_rate <- mean(diff(pop_sp$pop_est) / head(pop_sp$pop_est, -1))

# add projected rows for 2026 onward
forecast_years <- (last_year + 1):(last_year + ceiling(horizon / 12) + 1)
pop_projected <- data.frame(
  state    = "SP",
  year     = forecast_years,
  pop_est  = last_pop * (1 + growth_rate)^(forecast_years - last_year),
  observed = FALSE
)

pop_sp <- rbind(pop_sp, pop_projected)

# match to disease time vector and scale
population <- pop_sp$pop_est[match(as.integer(format(time, "%Y")), pop_sp$year)]
pop_scaled <- population / 1e6

#### step 4b: load stage2 forecast met (for step 7) ####
ds2 <- arrow::open_dataset(stage2_dir)
refs <- ds2 |>
  dplyr::distinct(.data$reference_datetime) |>
  dplyr::collect() |>
  dplyr::mutate(rd = as.Date(.data$reference_datetime))
fd <- as.Date(forecast_date)
ref_str <- format(if (fd %in% refs$rd) fd else max(refs$rd[refs$rd <= fd]))
reference_date <- as.Date(ref_str)
s2 <- ds2 |>
  dplyr::filter(.data$reference_datetime == ref_str) |>
  dplyr::collect()
met_stage2_monthly <- monthly_stage2(s2, sp_sites)

# Filter data
hasData <- !is.na(tMax) & !is.na(tMin) & !is.na(precip)
y       <- y[hasData]
tMax    <- tMax[hasData]
tMin    <- tMin[hasData]
precip  <- precip[hasData]
time    <- time[hasData]
pop_scaled <- pop_scaled[hasData]
n       <- length(y)

#### step 5: fit model to current data ####

# data list for JAGS
jags_data <- list(
  y      = as.integer(y),
  n      = n,
  tMax   = tMax,
  tMin   = tMin,
  precip = precip,
  pop    = pop_scaled,
  
  x_ic   = log(mean(y) + 1),
  tau_ic = 1,
  a_add  = 1,
  r_add  = 1,
  a_r    = 1,
  r_r    = 1
)

# JAGS model string
StateSpace = "
model{
  # Data Model (Negative Binomial)
  for(t in 1:n){
    y[t] ~ dnegbin(p[t], r)
    p[t] <- r / (r + mu[t])
    log(mu[t]) <- x[t]
  }
  # Process Model (AR(1) with climate covariates)
  for(t in 2:n){
    x[t] ~ dnorm(mu_proc[t], tau_add)
    mu_proc[t] <- alpha + phi * x[t-1] + 
                  beta_tMax * tMax[t] + 
                  beta_tMin * tMin[t] + 
                  beta_precip * precip[t] +
                  beta_pop * pop[t]
  }
  # Priors
  x[1] ~ dnorm(x_ic, tau_ic)
  tau_add ~ dgamma(a_add, r_add)
  r ~ dgamma(1, 0.01) # tightened from dgamma(0.001, 0.001) to improve convergence 
  alpha ~ dnorm(0, 0.01)
  phi ~ dnorm(0, 0.01)
  beta_tMax ~ dnorm(0, 0.01)
  beta_tMin ~ dnorm(0, 0.01)
  beta_precip ~ dnorm(0, 0.01)
  beta_pop ~ dnorm(0, 0.01)}"

# Define initial model state
nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y, length(y), replace=TRUE)
  init[[i]] <- list(
    tau_add = 1 / var(diff(log(y.samp + 1))),  
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
  variable.names = c("x",
                     "alpha", "phi",
                     "beta_tMax", "beta_tMin",
                     "beta_precip", "beta_pop",
                     "tau_add", "r"),
  n.iter = 10000)

# extract posterior matrix
out <- as.matrix(jags.out)


#### step 6: get initial conditions from last fitted state ####

# grab all columns that correspond to latent state x
x.cols <- grep("^x", colnames(out))

# notes!
# x[T] = last latent state = initial conditons for the forecast
# this is a distribution (one value per MCMC sample), not a single number !
x_T_all <- out[, max(x.cols)]

# subsample to n_ensemble members to propagate IC uncertainty into the forecast
idx    <- sample(nrow(out), n_ensemble, replace = FALSE)
x_T    <- x_T_all[idx]       # IC vector, length n_ensemble
params <- out[idx, ]          # matching parameter samples for each ensemble member

#### step 7: run ensemble forecast forward ####

# future months to forecast
future_months <- seq(
  from = max(time) %m+% months(1),
  by = "month",
  length.out = horizon)


# get future met drivers
# stage2 (31 members) × horizon
met_param <- (seq_len(n_ensemble) - 1L) %% 31L
mo <- as.Date(format(future_months,"%Y-%m-01"))

tMax_future   <- matrix(mean(tail(tMax, 12)),   nrow = horizon, ncol = n_ensemble)
tMin_future   <- matrix(mean(tail(tMin, 12)),   nrow = horizon, ncol = n_ensemble)
precip_future <- matrix(mean(tail(precip, 12)), nrow = horizon, ncol = n_ensemble)

# for (e in seq_len(n_ensemble)) {
#   sub <- dplyr::filter(met_stage2_monthly, .data$parameter == met_param[e])
#   m <- match(mo, sub$month)
#   tMax_future[, e]   <- sub$TMAX[m]
#   tMin_future[, e]   <- sub$TMIN[m]
#   precip_future[, e] <- sub$APCP[m]
# }
# get future population
pop_future <- pop_sp$pop_est[match(as.integer(format(future_months, "%Y")), pop_sp$year)] / 1e6

# initialize forecast storage
x_fc <- matrix(NA, nrow = horizon, ncol = n_ensemble)  # latent state
y_fc <- matrix(NA, nrow = horizon, ncol = n_ensemble)  # predicted cases

# run forecast for each ensemble member
for(e in 1:n_ensemble){
  
  x_prev <- x_T[e]        # ← initialize x_prev to the IC for this member
  x_fc[1, e] <- x_prev    # store it in the matrix too
  
  for(t in 2:horizon){
    mu_proc <- params[e, "alpha"] +
      params[e, "phi"] * x_prev +          # ← now correctly uses previous state
      params[e, "beta_tMax"]   * tMax_future[t, e] +
      params[e, "beta_tMin"]   * tMin_future[t, e] +
      params[e, "beta_precip"] * precip_future[t, e] +
      params[e, "beta_pop"]    * pop_future[t]
    
    x_fc[t, e] <- rnorm(1, mu_proc, sd = 1/sqrt(params[e, "tau_add"]))
    x_prev <- x_fc[t, e]   # ← carry state forward
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
fitted_ci  <- apply(mu_fitted, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
forecast_ci <- apply(y_fc, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

# build plot data frames
fitted_df <- tibble(
  month  = time,
  lo     = fitted_ci["2.5%",],
  median = fitted_ci["50%",],
  hi     = fitted_ci["97.5%",],
  obs    = y
)

forecast_df <- tibble(
  month  = future_months,
  lo     = forecast_ci["2.5%",],
  median = forecast_ci["50%",],
  hi     = forecast_ci["97.5%",]
)

# plot
# zoom to last year of historical + forecast
last_year_start <- max(time) %m-% months(12)

fitted_zoom <- filter(fitted_df, month >= last_year_start)

ggplot() +
  geom_ribbon(data = fitted_zoom, aes(x = month, ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.3) +
  geom_line(data = fitted_zoom, aes(x = month, y = median), color = "steelblue", linewidth = 0.8) +
  geom_point(data = fitted_zoom, aes(x = month, y = obs), color = "black", size = 1.5) +
  geom_ribbon(data = forecast_df, aes(x = month, ymin = lo, ymax = hi), fill = "tomato", alpha = 0.3) +
  geom_line(data = forecast_df, aes(x = month, y = median), color = "tomato", linewidth = 1) +
  geom_vline(xintercept = as.numeric(max(time)), linetype = "dashed", color = "grey40") +
  labs(title = "VL Forecast — São Paulo (last 12 months + forecast)", x = "Month", y = "Cases") +
  theme_bw()

#### step 9: write forecast file and submit ####

# build forecast dataframe
forecast_long <- tibble(
  model_id           = model_id,
  reference_datetime = reference_date,
  datetime           = rep(future_months, times = n_ensemble),
  site_id            = "SP",
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
s3_submit <- arrow::s3_bucket("bu4cast-ci-write/challenges/project_id=bu4cast/forecasts/", endpoint_override = s3_endpoint, anonymous = FALSE)  # need credentials to write

# upload file
#aws.s3::put_object(file = forecast_file, object = forecast_file, bucket = "bu4cast-ci-write",region = "", base_url = s3_endpoint)
