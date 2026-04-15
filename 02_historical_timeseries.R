# "02_historical_timeseries.R
# Script for fitting historical data using Bayesian state-space model
# Corresponds to 3/20/2026 "Pulling and Visualizing Data" Project Milestone


#### load packages and data ####
library(rjags)
library(ecoforecastR)
library(dplyr)

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
n <- nrow(model_site_monthly)
tMax <- model_site_monthly$tmax_c
tMin <- model_site_monthly$tmin_c
precip <- model_site_monthly$prec_mm
pop_scaled <- model_site_monthly$pop_scaled #scaled per 100k people

## remove the time points where driver data is unavailable
# find rows that DO have driver data
hasData <- !is.na(tMax) & !is.na(tMin) & !is.na(precip) & !is.na(pop_scaled)

# subset data and covariates (must keep lengths aligned and update n)
y <- y[hasData]
tMax <- tMax[hasData]
tMin <- tMin[hasData]
precip <- precip[hasData]
pop_scaled <- pop_scaled[hasData]
time <- time[hasData]
n <- length(y)

# Recreate data list for JAGS
jags_data <- list(
  y = as.integer(y),
  n = n,
  tMax = tMax,
  tMin = tMin,
  precip = precip,
  pop = pop_scaled,
  x_ic = log(mean(y, na.rm = TRUE) + 1),
  tau_ic = 1,
  
  a_add = 1,
  r_add = 1
)


#### Define models & priors as string to be passed into JAGS ####

StateSpace = "
model{
  
 #### Data Model (Negative Binomial)
  for(t in 1:n){
    y[t] ~ dnegbin(p[t], r)
    p[t] <- r / (r + mu[t])
    log(mu[t]) <- x[t]
  }
  
 #### Process Model (AR(1) with climate covariates)
for(t in 2:n){
  x[t] ~ dnorm(mu_proc[t], tau_add)
  
  mu_proc[t] <- alpha + phi * x[t-1] + 
  beta_tMax * tMax[t] + 
  beta_tMin * tMin[t] + 
  beta_precip * precip[t] +
  beta_pop * pop[t]
  }

  #### Priors
  x[1] ~ dnorm(x_ic, tau_ic) # prior for number of cases in first month
  
  tau_add ~ dgamma(a_add, r_add) # 
  
  # Dispersion parameter (size) for negative binomial
  r ~ dgamma(0.001, 0.001)
  
  
  alpha ~ dnorm(0, 0.01) # AR(1) intercept
  phi ~ dnorm(0, 0.01) # AR(1) coefficient = persistence 

  beta_tMax ~ dnorm(0, 0.01)  #temp max covariate
  beta_tMin ~ dnorm(0, 0.01)  #temp min covariate
  beta_precip ~ dnorm(0, 0.01) #precip. covariate
  beta_pop ~ dnorm(0, 0.01) #population covariate (NOT a prior for total population)
}
"
# Notes
# JAGS parameterizes NB with r (size/dispersion) and p (success probability)
# Ensure mu[t] > 0
# Convert mean -> probability


# Define initial model state
nchain = 3 # num. of MCMC chains
init <- list()

for(i in 1:nchain){
  y.samp = sample(y, length(y), replace=TRUE)
  
  init[[i]] <- list(
    tau_add = 1 / var(diff(log(y.samp + 1))),  
    r = 10   # starting guess for dispersion
  )
}


# Send info to JAGS
j.model <- jags.model(
  file = textConnection(StateSpace),
  data = jags_data,
  inits = init,
  n.chains = nchain
)

# MCMC chain samples
## burn-in
jags.out <- coda.samples(
  model = j.model,
  variable.names = c("tau_add","r"),
  n.iter = 2000
)

#diagnostics
plot(jags.out)
dic.samples(j.model, 2000)

# Larger samples
jags.out <- coda.samples(
  model = j.model,
  variable.names = c("x","tau_add","r"),
  n.iter = 10000
)

# Visualization
time.rng = c(1,length(time))
out <- as.matrix(jags.out)

x.cols <- grep("^x\\[", colnames(out))

ci <- apply(exp(out[, x.cols]), 2, quantile, c(0.025, 0.5, 0.975))

plot(time, ci[2,], type='n',
     ylim = range(y, na.rm = TRUE),
     ylab = "Cases",
     xlab = "Date",
     xlim = range(time))

date_span_days <- as.numeric(diff(range(time)))
if (date_span_days < 365 * 3) {
  axis.Date(1,
            at = seq(min(time), max(time), by = 'month'),
            format = "%Y-%m")
} else {
  axis.Date(1,
            at = seq(min(time), max(time), by = 'year'),
            format = "%Y")
}

ecoforecastR::ciEnvelope(
  time, ci[1,], ci[3,],
  col=ecoforecastR::col.alpha("lightBlue",0.75)
)

points(time, y, pch="+", cex=0.5)


# save MCMC output locally 
# saveRDS(jags.out, file = "mcmc_output.rds")
