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

# was going to use stage3 data (monthly, 30 ensembles) but it only has data from the past two months right now, so I have used the same daymet data that we used in the previous script as a placeholder right now...
# heres the code I wouldve used to bring that data in tho 
# s3 <- arrow::s3_bucket(paste0(drivers_path, "stage3/"), endpoint_override = s3_endpoint, anonymous = TRUE)
# met_hist <- arrow::open_dataset(s3) |> filter(site_id %in% sp_sites) |> collect()

daymet0 <- download_daymet(site = "SaoPaulo", start = 2007, end = 2025, internal = TRUE)
daymet <- daymet0$data

# make date col
daymet$date <- as.Date(paste(daymet$year, daymet$yday, sep = "-"), "%Y-%j")
daymet$month <- as.Date(format(daymet$date, "%Y-%m-01"))

# aggregate to monthly
tMin_monthly   <- aggregate(daymet$tmin..deg.c., by = list(daymet$month), mean, na.rm = TRUE)
tMax_monthly   <- aggregate(daymet$tmax..deg.c., by = list(daymet$month), mean, na.rm = TRUE)
precip_monthly <- aggregate(daymet$prcp..mm.day., by = list(daymet$month), sum, na.rm = TRUE)

time_month <- as.Date(format(time, "%Y-%m-01"))

# extract and align to disease time vector
tMin   <- tMin_monthly[,2][match(time_month, tMin_monthly[,1])]
tMax   <- tMax_monthly[,2][match(time_month, tMax_monthly[,1])]
precip <- precip_monthly[,2][match(time_month, precip_monthly[,1])]

# remove time points where climate data is unavailable
hasData <- !is.na(tMax) & !is.na(tMin) & !is.na(precip)
y      <- y[hasData]
tMax   <- tMax[hasData]
tMin   <- tMin[hasData]
precip <- precip[hasData]
time   <- time[hasData]
n      <- length(y)

met_hist <- arrow::open_dataset(s3) |> filter(site_id %in% sp_sites) |> collect()

# Align met drivers to disease time vector, handle NAs


#### step 3: load future met drivers (stage2) ####

stage2 <- arrow::s3_bucket(paste0(drivers_path, "stage2/"),endpoint_override = s3_endpoint,anonymous = TRUE)
met_future <- arrow::open_dataset(s3_stage2) |> filter(site_id %in% sp_sites) |> collect()

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

#### step 5: fit model to current data ####

# data list for JAGS
data <- list(
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
  data = data,
  inits = init,
  n.chains = nchain)

# burn-in
jags.out <- coda.samples(
  model = j.model,
  variable.names = c("tau_add", "r"),
  n.iter = 2000)

# diagnostics
plot(jags.out)
dic.samples(j.model, 2000)

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
  from = as.Date(format(Sys.Date(), "%Y-%m-01")),
  by = "month",
  length.out = horizon)


# get future met drivers
# current sampling from daymet as a placeholder until stage2 is available
tMax_future   <- sample(tMax, horizon, replace = TRUE)
tMin_future   <- sample(tMin, horizon, replace = TRUE)
precip_future <- sample(precip, horizon, replace = TRUE)

# get future population
pop_future <- pop_sp$pop_est[match(as.integer(format(future_months, "%Y")), pop_sp$year)] / 1e6

# initialize forecast storage
x_fc <- matrix(NA, nrow = horizon, ncol = n_ensemble)  # latent state
y_fc <- matrix(NA, nrow = horizon, ncol = n_ensemble)  # predicted cases

# run forecast for each ensemble member
for(e in 1:n_ensemble){
  
  x_fc[1, e] <- x_T[e]  # initialize from IC
  
  for(t in 2:horizon){
    # process model
    mu_proc <- params[e, "alpha"] + 
      params[e, "phi"] * x_fc[t-1, e] +
      params[e, "beta_tMax"] * tMax_future[t] +
      params[e, "beta_tMin"] * tMin_future[t] +
      params[e, "beta_precip"] * precip_future[t] +
      params[e, "beta_pop"] * pop_future[t]
    
    # add process error
    x_fc[t, e] <- rnorm(1, mu_proc, sd = 1/sqrt(params[e, "tau_add"]))
  }
  
  # observation model -- convert latent state to predicted cases
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
# reference_datetime = forecast_date (today's run date)
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



#### step 9: write forecast file and submit ####

# build forecast dataframe
forecast_long <- tibble(
  model_id           = model_id,
  reference_datetime = forecast_date,
  datetime           = rep(future_months, times = n_ensemble),
  site_id            = "SP",
  family             = "ensemble",
  parameter          = rep(1:n_ensemble, each = horizon),
  variable           = "cases",
  prediction         = as.vector(y_fc)
)

head(forecast_long)

# filename convention: "<category>-<reference_datetime>-<model_id>.csv.gz"
forecast_file <- paste0("disease-", forecast_date, "-", model_id, ".csv.gz")

# write compressed csv
write_csv(forecast_long, forecast_file)
file.exists(forecast_file) # did it work?

# submit to BU challenge S3 bucket 
# COMMENTED OUT UNTIL WE ARE READY!
# s3_submit <- arrow::s3_bucket("bu4cast-ci-write/challenges/project_id=bu4cast/forecasts/", endpoint_override = s3_endpoint, anonymous = FALSE  # need credentials to write

# upload file
#aws.s3::put_object(file = forecast_file, object = forecast_file, bucket = "bu4cast-ci-write",region = "", base_url = s3_endpoint)
