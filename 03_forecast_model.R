# 03_forecast_model.R
# Iterative ensemble forecast for tropical disease (VL cases, all 8 sites)
# Corresponds to 4/20/2026 "Data Assimilation" Milestone

#### outline ####
#   1. load target data (all 8 sites, 2007-2025)
#   2. outer DA loop: Sep 2024 - Jul 2025, one month at a time
#       a. filter to data available up to that month
#       b. for each site: fit state-space model on all data up to today
#       c. get last latent state from the fit to use as IC for next steps
#       d. forecast forward horizon months
#       e. format in EFI standard format
#   3. combine all sites/months, write forecast file, submit
#   4. visualize fitted + forecast for all 8 sites in one panel

# NOTE on "starts from previous forecast":
#   This does NOT mean loading yesterday's output file.
#   It means each daily run re-fits the model with the latest data,
#   and uses X[T] (the last state of that fit) as the IC for the forecast.
#   Since our model is fast enough, we re-fit from scratch each run.


#### step 0: load packages and settings ####
library(rjags)
library(ecoforecastR)
library(tidyverse)
library(lubridate)
library(arrow)
library(patchwork)

source("00_clean_and_plot_monthly_inputs.R")

model_id      <- "BU_visceral_leishmaniacs"
forecast_date <- Sys.Date()
n_ensemble    <- 100
horizon       <- 2  # months forward 

s3_endpoint  <- "minio-s3.apps.shift.nerc.mghpcc.org"
target_url   <- "https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-read/challenges/project_id=bu4cast/targets/tropical-disease-targets.csv"
drivers_path <- "bu4cast-ci-read/challenges/project_id=bu4cast/drivers/"

# our 8 target sites (one per climate region)
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

# data asimilation holdout period: one step per month
assimilation_months <- seq(as.Date("2024-08-01"), as.Date("2025-07-01"), by = "month")

# helper function to get monthly forecast met from stage2
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

# JAGS model string
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

#### step 1: load data for all sites####
# load full historical data (2007-2025) 
full_df <- monthly_data(site_ids   = site_ids, start_date = as.Date("2007-01-01"), end_date   = as.Date("2025-08-01"))$monthly_merged

#### step 2: load stage2 met data ####
# changed to load stage2 once outside the loop 
stage2         <- arrow::s3_bucket(paste0(drivers_path, "stage2/"), endpoint_override = s3_endpoint, anonymous = TRUE)

# subset before loading to fix memory issue
met_future_raw <- arrow::open_dataset(stage2) |> 
  dplyr::filter(site_id %in% site_ids) |>
  dplyr::collect()

#### step 3: data asimillation loop ####
# this is where i added my big outer loop! each iteration adds one more month of observed data

all_forecasts <- list()
plot_list     <- list()

for (end_date in as.list(assimilation_months)) {
  end_date <- as.Date(end_date)
  
  # filter to data available up to this DA step
  df <- full_df %>% filter(month <= end_date)
  
  # run model independently for each site
  for (site in site_ids) {
    
    model_site_id      <- site  
    model_site_monthly <- df %>%
      filter(site_id == site) %>%
      arrange(month)
    
    if (nrow(model_site_monthly) < 12) next
    
    y          <- model_site_monthly$observation
    time       <- model_site_monthly$month
    temp       <- (model_site_monthly$tmin_c + model_site_monthly$tmax_c) / 2
    precip     <- model_site_monthly$prec_mm
    pop_scaled <- model_site_monthly$pop_scaled
    
    # remove months with missing drivers (e.g. 2023 population gap)
    hasData    <- !is.na(temp) & !is.na(precip) & !is.na(pop_scaled)
    y          <- y[hasData]
    time       <- time[hasData]
    temp       <- temp[hasData]
    precip     <- precip[hasData]
    pop_scaled <- pop_scaled[hasData]
    n_obs      <- length(y)  # reset AFTER filtering
    
    if (n_obs < 12) next
    
    forecast_start <- max(time) %m+% months(1)
    future_months  <- seq(forecast_start, by = "month", length.out = horizon)
    
    #### stage2 met for this site ####
    # i changed this to use pre-loaded met_future_raw instead of re-downloading
    met_future <- met_future_raw
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
    
    #### population projection ####
    pop_site <- full_df %>%
      filter(site_id == site) %>%
      distinct(year = as.integer(format(month, "%Y")), pop_est = pop_scaled * 1e6) %>%
      arrange(year)
    
    last_year   <- max(pop_site$year)
    last_pop    <- pop_site$pop_est[pop_site$year == last_year][1]
    growth_rate <- mean(diff(pop_site$pop_est) / head(pop_site$pop_est, -1))
    
    forecast_years <- (last_year + 1):(last_year + ceiling(horizon / 12) + 1)
    pop_projected  <- data.frame(
      year    = forecast_years,
      pop_est = last_pop * (1 + growth_rate)^(forecast_years - last_year)
    )
    pop_site <- rbind(pop_site, pop_projected)
    
    #### fit JAGS ####
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
                             c("x", "alpha", "phi", "beta_temp", "beta_precip", "beta_pop", "tau_add", "r"),
                             n.iter = 10000)
    out <- as.matrix(jags.out)
    
    #### initial conditions from last fitted state ####
    x.cols <- grep("^x\\[", colnames(out))
    idx    <- sample(nrow(out), n_ensemble, replace = nrow(out) < n_ensemble)
    x_T    <- out[idx, max(x.cols)]
    params <- out[idx, ]
    
    #### future met drivers ####
    met_monthly_forecast <- monthly_met_future(met_future, model_site_id)
    n_met_members        <- dplyr::n_distinct(met_monthly_forecast$parameter)
    met_param            <- (seq_len(n_ensemble) - 1) %% max(n_met_members, 1L)
    mo                   <- as.Date(format(future_months, "%Y-%m-01"))
    mo_df                <- tibble(month = mo)
    
    temp_future   <- matrix(NA_real_, nrow = horizon, ncol = n_ensemble)
    precip_future <- matrix(NA_real_, nrow = horizon, ncol = n_ensemble)
    
    if (n_met_members == 0L) {
      temp_future[]   <- mean(temp,   na.rm = TRUE)
      precip_future[] <- mean(precip, na.rm = TRUE)
    } else {
      for (e in seq_len(n_ensemble)) {
        sub <- dplyr::filter(met_monthly_forecast, parameter == met_param[e]) |>
          dplyr::mutate(month = as.Date(as.character(month)))
        z                  <- dplyr::left_join(mo_df, sub, by = "month")
        temp_future[, e]   <- z$TMP
        precip_future[, e] <- z$APCP
      }
      clim_temp   <- mean(temp,   na.rm = TRUE)
      clim_precip <- mean(precip, na.rm = TRUE)
      for (e in seq_len(n_ensemble)) {
        temp_future[is.na(temp_future[, e]),     e] <- clim_temp
        precip_future[is.na(precip_future[, e]), e] <- clim_precip
      }
    }
    
    pop_future <- pop_site$pop_est[match(as.integer(format(future_months, "%Y")), pop_site$year)] / 1e6
    pop_future[is.na(pop_future)] <- tail(pop_scaled[!is.na(pop_scaled)], 1)
    #### forecast ####
    x_fc <- y_fc <- matrix(NA, nrow = horizon, ncol = n_ensemble)
    
    for (e in seq_len(n_ensemble)) {
      x_prev <- x_T[e]
      for (t in seq_len(horizon)) {
        mu_proc    <- params[e, "alpha"] +
          params[e, "phi"]         * x_prev           +
          params[e, "beta_temp"]   * temp_future[t, e]   +
          params[e, "beta_precip"] * precip_future[t, e] +
          params[e, "beta_pop"]    * pop_future[t]
        x_fc[t, e] <- rnorm(1, mu_proc, 1 / sqrt(params[e, "tau_add"]))
        x_prev     <- x_fc[t, e]
      }
      for (t in seq_len(horizon)) {
        mu_t       <- exp(x_fc[t, e])
        p_t        <- params[e, "r"] / (params[e, "r"] + mu_t)
        y_fc[t, e] <- rnbinom(1, size = params[e, "r"], prob = p_t)
      }
    }
    
    #### store EFI output ####
    # make sure in proper EFi format....
    all_forecasts[[paste0(site, "_", end_date)]] <- tibble(
      model_id           = model_id,
      reference_datetime = end_date,
      datetime           = rep(future_months, times = n_ensemble),
      site_id            = site,
      family             = "ensemble",
      parameter          = rep(1:n_ensemble, each = horizon),
      variable           = "cases",
      prediction         = as.vector(y_fc)
    )
    
    #### visualize !! ####
    # this code only builds plots on the last DA step because if not R crashes...
    if (end_date == max(assimilation_months)) {
      x_fitted_cols <- matrix(out[idx, grep("^x\\[", colnames(out))], nrow = n_ensemble)
      fitted_ci     <- apply(exp(x_fitted_cols), 2, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
      forecast_ci_lo  <- apply(y_fc, 1, function(v) quantile(v, 0.025, na.rm = TRUE))
      forecast_ci_med <- apply(y_fc, 1, function(v) quantile(v, 0.5,   na.rm = TRUE))
      forecast_ci_hi  <- apply(y_fc, 1, function(v) quantile(v, 0.975, na.rm = TRUE))
      
      fitted_df <- tibble(
        month  = time,
        lo     = fitted_ci["2.5%",],
        median = fitted_ci["50%",],
        hi     = fitted_ci["97.5%",],
        obs    = y
      )
      forecast_df <- tibble(
        month  = future_months,
        lo     = forecast_ci_lo,
        median = forecast_ci_med,
        hi     = forecast_ci_hi
      )
      
      last_year_start <- max(time) %m-% months(12)
      fitted_zoom     <- filter(fitted_df, month >= last_year_start)
      x_rng           <- range(c(fitted_zoom$month, forecast_df$month), na.rm = TRUE)
      
      # plot aesthetics
      plot_list[[as.character(site)]] <- ggplot() +
        geom_ribbon(data = fitted_zoom,
                    aes(x = month, ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.3) +
        geom_line(data = fitted_zoom,
                  aes(x = month, y = median), color = "steelblue", linewidth = 0.8) +
        geom_point(data = fitted_zoom,
                   aes(x = month, y = obs), color = "black", size = 1.5) +
        geom_ribbon(data = forecast_df,
                    aes(x = month, ymin = lo, ymax = hi), fill = "tomato", alpha = 0.3) +
        geom_line(data = forecast_df,
                  aes(x = month, y = median), color = "tomato", linewidth = 1) +
        geom_vline(xintercept = as.numeric(max(time)), linetype = "dashed", color = "grey40") +
        scale_x_date(limits = x_rng, date_breaks = "2 months", date_labels = "%b %Y") +
        labs(title = site_names[as.character(site)], x = "Month", y = "Cases") +
        theme_bw()
    }
    
    # added this to clear memory after each site to prevent R from crashing
    rm(jags.out, out, j.model, x_fc, y_fc, params)
    gc()
    
  }  # end site loop
}  # end DA loop

#### step 4: combine and write forecast ####
combined_forecast <- bind_rows(all_forecasts)
forecast_file     <- paste0("disease-", Sys.Date(), "-", model_id, ".csv.gz")
write_csv(combined_forecast, forecast_file)

#### step 5: visualize ####
# combine all 8 site plots into one panel and save
dir.create("figures", showWarnings = FALSE)

panel <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(
    title    = "VL Forecast — All Sites",
    subtitle = paste0("Final DA step: ", max(assimilation_months))
  )

ggsave("figures/forecast_all_sites.png", plot = panel, width = 14, height = 20, dpi = 300)
print(panel)

#### step 6: submit ####
aws.s3::put_object(file     = forecast_file, object   = forecast_file, bucket   = "bu4cast-ci-write",region   = "",base_url = s3_endpoint)
