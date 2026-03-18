# "01_daily_download_and_vis.R'
# Script for daily target data download and visualization
# Corresponds to 2/20/2026 "Pulling and Visualizing Data" Project Milestone

#### part 1: load packages ####
# install.packages("readr")
# install.packages("ggplot2")
library(readr)
library(ggplot2)

#### part 2: bring in data ####
disease_url = 'https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-read/challenges/project_id=bu4cast/targets/tropical-disease-targets.csv'

disease_targets <- read.csv(disease_url)

#### part 3: plot time-series visualization of data ####

# split observations by time to create timeseries for each site
obs_by_time <- split(disease_targets$observation, disease_targets$datetime) 

# compute monthly quantiles
# code borrowed from EF_Activities/Exercise_02_Logistic
n.stats <- sapply(obs_by_time, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

time <- as.Date(names(obs_by_time)) # convert from character to Date for plotting
ylo  <- n.stats[1, ]
ymed <- n.stats[2, ]
yhi  <- n.stats[3, ]

# function to draw a confidence interval ribbon between ylo and yhi across time points in x
ciEnvelope <- function(x, ylo, yhi, col = rgb(0,0,0,0.2), ...) {
  polygon(c(x, rev(x), x[1]),
          c(ylo, rev(yhi), ylo[1]),
          border = NA, col = col, ...)
}

par(mfrow = c(3, 1),mar = c(3, 4, 3, 1))

# to see raw points and get a sense for distribution across sites
plot(as.Date(disease_targets$datetime), disease_targets$observation,
     ylab = "VL Cases", col=rgb(0,0,0,0.2), pch = 16,
     main = "Raw monthly observations")

# to get a sense for variability across sites (zoom on median+CI)
plot(time, ymed, type = "l",
     ylim = range(c(ylo, yhi)),
     ylab = "VL Cases",
     main = "Median trajectory across sites with 95% CI")

ciEnvelope(time, ylo, yhi) # add confidence interval
lines(time, ymed, lwd = 2) # redraw median on top

# total cases per monthly across sites (more standard for this kind of data!)
total_cases <- sapply(obs_by_time, sum)

barplot(total_cases, space=0, names.arg=format(time, "%Y   "), las=1,
        ylab="VL cases",
        main="Montly total VL cases across sites")