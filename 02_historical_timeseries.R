# "02_historical_timeseries.R
# Script for fitting historical data using Bayesian state-space model
# Corresponds to 3/20/2026 "Pulling and Visualizing Data" Project Milestone



#### load packages and data ####
library(rjags)
library(daymetr)
library(ecoforecastR)

# get disease datat
disease_url = 'https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-read/challenges/project_id=bu4cast/targets/tropical-disease-targets.csv'

disease_targets = read.csv(disease_url)

# subset data to only Sao Paolo sites
sp <- subset(disease_targets, site_id >= 350000 & site_id < 360000)  # SP = 35xxxx IBGE prefix
sp_monthly <- aggregate(observation ~ datetime, data = sp, sum, na.rm = TRUE) #aggregate to monthly

# Format input data
time = as.Date(sp_monthly$datetime)
y = sp_monthly$observation

plot(time, y, type='l', ylab="Cases", lwd=2)

# get weather data for Sao Paulo (we may need to double check this site to make sure it actually overlaps with ours)
# we still dont have our actual met data, so this is our current best bet...
daymet0 <- download_daymet(site = "SaoPaulo", start = 2007, end = 2025, internal = TRUE)
daymet <- daymet0$data

# make date col
daymet$date <- as.Date(paste(daymet$year, daymet$yday, sep = "-"), "%Y-%j")

# extract temp and precip data
tMin <- daymet$tmin..deg.c.[match(time, daymet$date)]
tMax <- daymet$tmax..deg.c.[match(time, daymet$date)]
precip <- daymet$prcp..mm.day.[match(time, daymet$date)]

## remove the time points where climate data is unavailable
# find rows that DO have climate data
hasData <- !is.na(tMax) & !is.na(tMin) & !is.na(precip)

# subset data and covariates
y <- y[hasData]
tMax <- tMax[hasData]
tMin <- tMin[hasData]
precip <- precip[hasData]
time <- time[hasData]

# Update length for JAGS
n <- length(y)

## add population data! 
pop_path <- file.path("data", "state_pop_by_year.csv")
pop_df <- read.csv(pop_path)
pop_sp <- subset(pop_df, state == "SP")
pop_sp <- pop_sp[order(pop_sp$year), ]

# repeated monthly within year (to match disease monthly cadence)
population <- pop_sp$pop_est[match(as.integer(format(time, "%Y")), pop_sp$year)]

# divided by 1 mil so beta_pop is per mil people (similar scale to climate covariates)
pop_scaled <- population / 1e6

# Recreate data list for JAGS
data <- list(
  y = as.integer(y),
  n = n,
  tMax = tMax,
  tMin = tMin,
  precip = precip,
  pop = pop_scaled,
  x_ic = log(mean(y) + 1),
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
  data = data,
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
  variable.names = c("x","tau_add","r"),  #let's look at the other coefficients too at some point
  n.iter = 10000
)

# Visualization (needs fixed now)
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
