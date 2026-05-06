# 05_model_assessment.R
# Model assessment: comparing our null climatological vs dynamic state-space model using Continuous Ranked Probability Score (CRPS), following Activity 11


# Notes
# CRPS = (1/m) * sum|Xi - y| - (1/2m^2) * sum_i sum_j |Xi - Xj|
# where y = observed, Xi = ensemble members
# Lower CRPS = better forecast

#### load packages ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scoringRules) 
library(lubridate)


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

assimilation_months <- seq(as.Date("2024-08-01"), as.Date("2025-07-01"), by = "month")
horizon <- 2

dir.create("figures", showWarnings = FALSE)


#### load forecasts ####

# --- null model ---
null_fc <- readRDS("null_forecast_ensemble.rds")

# --- dynamic model ---
dyn_files <- list.files(".", pattern = "^disease-.*BU_visceral_leishmaniacs\\.csv\\.gz$",
                        full.names = TRUE)
if (length(dyn_files) == 0) {
  stop("No dynamic forecast file found. Run 03_forecast_model.R first.")
}
dyn_file <- sort(dyn_files, decreasing = TRUE)[1]

dyn_fc <- readr::read_csv(dyn_file, show_col_types = FALSE) %>%
  mutate(datetime = as.Date(datetime),
         reference_datetime = as.Date(reference_datetime))


#### load observations ####
obs_raw <- read.csv(
  "https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-read/challenges/project_id=bu4cast/targets/tropical-disease-targets.csv"
) %>%
  mutate(datetime = as.Date(substr(datetime, 1, 10))) %>%
  filter(site_id %in% site_ids) %>%
  select(site_id, datetime, observation)


#### helper function to compute CRPS ####
compute_crps <- function(forecast_df, obs_df, model_label) {
  forecast_df %>%
    inner_join(obs_df, by = c("site_id", "datetime")) %>%
    filter(!is.na(observation)) %>%
    group_by(site_id, datetime, reference_datetime, observation) %>%
    summarise(
      crps_val = scoringRules::crps_sample(
        y   = unique(observation),
        dat = prediction
      ),
      n_members = n(),
      .groups = "drop"
    ) %>%
    mutate(
      model      = model_label,
      horizon_mo = as.integer(
        round(as.numeric(difftime(datetime, reference_datetime, units = "days")) / 30.44)
      ),
      cal_month  = as.integer(format(datetime, "%m")),
      site_label = site_names[as.character(site_id)]
    )
}

crps_null <- compute_crps(null_fc, obs_raw, "Null")
crps_dyn  <- compute_crps(dyn_fc,  obs_raw, "Dynamic")
crps_all  <- bind_rows(crps_null, crps_dyn)


#### summary table ####
summary_tbl <- crps_all %>%
  group_by(model, site_label) %>%
  summarise(
    mean_CRPS   = round(mean(crps_val, na.rm = TRUE), 2),
    median_CRPS = round(median(crps_val, na.rm = TRUE), 2),
    n_forecasts = n(),
    .groups = "drop"
  ) %>%
  arrange(site_label, model)

print(summary_tbl, n = Inf)
readr::write_csv(summary_tbl, "figures/crps_summary_table.csv")


#### overall mean CRPS per model ####
overall <- crps_all %>%
  group_by(model) %>%
  summarise(mean_CRPS = round(mean(crps_val, na.rm = TRUE), 3), .groups = "drop")
print(overall)


#### skill score relative to null ####
skill_tbl <- summary_tbl %>%
  select(model, site_label, mean_CRPS) %>%
  pivot_wider(names_from = model, values_from = mean_CRPS) %>%
  rename(crps_null = `Null`,
         crps_dyn  = `Dynamic`) %>%
  mutate(skill = round(1 - crps_dyn / crps_null, 3))

print(skill_tbl)
readr::write_csv(skill_tbl, "figures/skill_scores.csv")

# standardize model labels
crps_all <- crps_all %>%
  mutate(model = recode(model,
                        "Null (climatology)"      = "Null",
                        "Dynamic (AR1 + climate)" = "Dynamic"))

summary_tbl <- summary_tbl %>%
  mutate(model = recode(model,
                        "Null (climatology)"      = "Null",
                        "Dynamic (AR1 + climate)" = "Dynamic"))

colors <- c("Null" = "#888888", "Dynamic" = "#2E75B6")


#### plot 1: CRPS over time ####
a1 <- crps_all %>%
  group_by(model, datetime) %>%
  summarise(mean_crps = mean(crps_val, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = datetime, y = mean_crps, color = model)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  scale_color_manual(values = colors) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(x = NULL, y = "Mean CRPS", color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 30, hjust = 1))


#### plot 2: skill score bar chart ####
a2 <- ggplot(skill_tbl, aes(x = reorder(site_label, skill), y = skill, fill = skill > 0)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  scale_fill_manual(values = c("TRUE" = "#2E75B6", "FALSE" = "#E74C3C"), guide = "none") +
  labs(x = NULL, y = "Skill Score") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#### plot 3: Bayesian p-values (quantile calibration) ####
# For each dynamic forecast, compute what quantile the observation falls in
# flat histogram = well calibrated
# U-shaped = overconfident (intervals too narrow)
# hump at 0.5 = underconfident (intervals too wide)

pval_dyn <- dyn_fc %>%
  inner_join(obs_raw, by = c("site_id", "datetime")) %>%
  filter(!is.na(observation)) %>%
  group_by(site_id, datetime) %>%
  summarise(
    pval = mean(observation > prediction),
    .groups = "drop"
  )

a3 <- ggplot(pval_dyn, aes(x = pval)) +
  geom_histogram(bins = 10, fill = "#2E75B6", color = "white", alpha = 0.85) +
  geom_hline(yintercept = nrow(pval_dyn) / 10,
             linetype = "dashed", color = "grey30") +
  labs(x = "Quantile of observation within ensemble", y = "Count") +
  theme_bw(base_size = 12)


#### plot 4: CRPS vs observed cases ####
a4 <- crps_all %>%
  ggplot(aes(x = observation, y = crps_val, color = model)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  scale_color_manual(values = colors) +
  labs(x = "Observed VL cases", y = "CRPS", color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")


#### combined 4-panel ####
plot_panel <- (a1 | a2) / (a3 | a4)

ggsave("figures/model_assessment.png", plot_panel,
       width = 14, height = 10, dpi = 300)
print(plot_panel)


#### 8-panel site skill figure ####
site_time_crps <- crps_all %>%
  group_by(model, site_label, datetime) %>%
  summarise(mean_crps = mean(crps_val, na.rm = TRUE), .groups = "drop")

skill_labels <- skill_tbl %>%
  mutate(skill_label = paste0("Skill score = ", skill))

site_order <- site_names[as.character(site_ids)]
site_time_crps <- site_time_crps %>%
  mutate(site_label = factor(site_label, levels = site_order))
skill_labels <- skill_labels %>%
  mutate(site_label = factor(site_label, levels = site_order))

site_plots <- lapply(site_order, function(sl) {
  skill_val   <- skill_labels %>% filter(site_label == sl) %>% pull(skill)
  skill_str   <- if (length(skill_val) == 0 || is.na(skill_val)) "N/A" else paste0("Skill = ", skill_val)
  skill_color <- if (!is.na(skill_val) && skill_val > 0) "#2E75B6" else "#E74C3C"
  
  site_time_crps %>%
    filter(site_label == sl) %>%
    ggplot(aes(x = datetime, y = mean_crps, color = model)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    scale_color_manual(values = colors) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b %y") +
    labs(title = sl, subtitle = skill_str, x = NULL, y = "Mean CRPS", color = NULL) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title      = element_text(size = 9,  face = "bold"),
      plot.subtitle   = element_text(size = 8,  color = skill_color, face = "bold"),
      axis.text.x     = element_text(angle = 30, hjust = 1, size = 7),
      axis.text.y     = element_text(size = 7),
      plot.margin     = margin(4, 6, 4, 6)
    )
})

legend_plot <- ggplot(
  data.frame(model = c("Null", "Dynamic"), x = 1, y = 1),
  aes(x = x, y = y, color = model)
) +
  geom_line() +
  scale_color_manual(values = colors, name = NULL) +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 11))

legend_grob <- cowplot::get_legend(legend_plot)

site_panel_grid <- wrap_plots(site_plots, ncol = 4) +
  plot_annotation(
    title    = "Forecast Skill by Site",
    subtitle = "skill = (1 − CRPS_dynamic / CRPS_null)",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

final_site_panel <- cowplot::plot_grid(
  site_panel_grid, legend_grob,
  ncol = 1, rel_heights = c(1, 0.04)
)

ggsave("figures/skill_by_site.png", final_site_panel,
       width = 16, height = 9, dpi = 300)