# 04_null_model.R
# Climatological null model for VL forecasting
# for each calendar month (1-12), predict cases using the historical mean and SD of observations at that site/month combination. this captures seasonal climatology without any process dynamics


#### load packages ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

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

n_ensemble   <- 100
# DA holdout window: same as dynamic model for fair comparison
train_end    <- as.Date("2024-08-01")
assimilation_months <- seq(as.Date("2024-08-01"), as.Date("2025-07-01"), by = "month")
horizon      <- 2


#### load data ####
d      <- monthly_data(site_ids = site_ids,
                       start_date = as.Date("2007-01-01"),
                       end_date   = as.Date("2025-08-01"))
df     <- d$monthly_merged %>%
  mutate(cal_month = as.integer(format(month, "%m")))  # 1-12 calendar month


#### fit climatology per site & calendar month ####
# Training data: everything up to the holdout window
clim_stats <- df %>%
  filter(month < train_end) %>%
  group_by(site_id, cal_month) %>%
  summarise(
    clim_mean = mean(observation, na.rm = TRUE),
    clim_sd   = sd(observation,   na.rm = TRUE),
    clim_n    = sum(!is.na(observation)),
    .groups   = "drop"
  ) %>%
  # if only 1 obs, sd is NA → fall back to Poisson-like sqrt(mean)
  mutate(clim_sd = ifelse(is.na(clim_sd) | clim_sd == 0,
                          sqrt(pmax(clim_mean, 1)),
                          clim_sd))


#### iterative null forecast (mirrors DA loop in 03_forecast_model.R) ####
# For each DA step and each site, draw n_ensemble samples from
# Normal(clim_mean, clim_sd) truncated at 0 for the forecast months.
# This gives a proper ensemble so we can compute CRPS fairly.

all_null_forecasts <- list()

for (end_date in as.list(assimilation_months)) {
  end_date <- as.Date(end_date)
  
  for (site in site_ids) {
    
    forecast_start <- end_date %m+% months(1)  # one month ahead of DA cutoff
    future_months  <- seq(forecast_start, by = "month", length.out = horizon)
    
    for (h in seq_along(future_months)) {
      fm         <- future_months[h]
      cal_mo     <- as.integer(format(fm, "%m"))
      
      stats <- clim_stats %>%
        filter(site_id == site, cal_month == cal_mo)
      
      if (nrow(stats) == 0) next
      
      # draw ensemble — truncate at 0 (cases can't be negative)
      draws <- pmax(0, round(
        rnorm(n_ensemble, mean = stats$clim_mean, sd = stats$clim_sd)
      ))
      
      all_null_forecasts[[paste0(site, "_", end_date, "_h", h)]] <- tibble(
        model_id           = "null_climatological",
        reference_datetime = end_date,
        datetime           = fm,
        site_id            = site,
        family             = "ensemble",
        parameter          = 1:n_ensemble,
        variable           = "cases",
        prediction         = draws
      )
    }
  }
}

null_forecast_df <- bind_rows(all_null_forecasts)

# save for use in model assessment script
saveRDS(null_forecast_df, file = "null_forecast_ensemble.rds")
message("Null forecast saved: null_forecast_ensemble.rds")


#### visualize null climatology ####
# Show seasonal mean ± 1 SD per site to confirm the null makes sense

clim_plot_df <- clim_stats %>%
  mutate(
    site_label = site_names[as.character(site_id)],
    lo = pmax(0, clim_mean - clim_sd),
    hi = clim_mean + clim_sd
  )

p_clim <- ggplot(clim_plot_df,
                 aes(x = cal_month, y = clim_mean,
                     color = site_label, fill = site_label,
                     group = site_label)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_x_continuous(breaks = 1:12,
                     labels = month.abb) +
  labs(x        = "Month",
    y        = "VL Cases (historical mean)",
    color    = "Site",
    fill     = "Site"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

print(p_clim)
dir.create("figures", showWarnings = FALSE)
ggsave("figures/null_climatology.png", p_clim, width = 12, height = 5, dpi = 300)
