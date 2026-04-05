# Time series of VL cases for Sao Paulo (SP)
# Dataset has no `state` column; SP is identified by IBGE prefix 35 in `site_id`

library(ggplot2)
library(readr)

disease_url <- "https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-read/challenges/project_id=bu4cast/targets/tropical-disease-targets.csv"
d <- read.csv(disease_url, stringsAsFactors = FALSE)
output_dir <- "/projectnb/dietzelab/skanee/ee585"

# Keep SP rows using 35xxxx site_id and convert timestamp -> date
d$date <- as.Date(sub(" 00:00:00.000000Z", "", d$datetime, fixed = TRUE))
sp <- d[d$site_id >= 350000 & d$site_id < 360000, ]
sp$site_id <- factor(sp$site_id)

if (nrow(sp) == 0) stop("No rows found for SP (site_id 350000-359999).")

# Total VL cases per month across all SP sites (statewide time series)
sp_monthly <- aggregate(observation ~ date, data = sp, sum, na.rm = TRUE)
sp_monthly$date <- as.Date(sp_monthly$date)
names(sp_monthly)[names(sp_monthly) == "observation"] <- "VL_cases"
sp_monthly <- sp_monthly[order(sp_monthly$date), c("date", "VL_cases")]

# Save statewide SP monthly time series for downstream modeling
write.csv(sp_monthly, file.path(output_dir, "VL_timeseries_SP_monthly.csv"), row.names = FALSE)

# ---- Plot 1: Total VL cases per month (all SP sites combined) ----
p_total <- ggplot(sp_monthly, aes(x = date, y = VL_cases)) +
  geom_line(linewidth = 0.8, color = "steelblue") +
  geom_point(size = 1.5, color = "steelblue", alpha = 0.8) +
  labs(
    x = "Date",
    y = "VL cases",
    title = "Total VL cases per month - SP (all site_ids)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(output_dir, "VL_timeseries_SP_total.png"), p_total, width = 8, height = 4, dpi = 150)
print(p_total)

# ---- Plot 2: One line per site_id (overlaid, semi-transparent) ----
p_sites <- ggplot(sp, aes(x = date, y = observation, group = site_id)) +
  geom_line(alpha = 0.3, linewidth = 0.4, color = "darkblue") +
  geom_line(
    data = sp_monthly,
    aes(x = date, y = VL_cases, group = 1),
    inherit.aes = FALSE,
    linewidth = 1,
    color = "red"
  ) +
  labs(
    x = "Date",
    y = "VL cases",
    title = "VL cases by site_id in SP (red = statewide total)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(output_dir, "VL_timeseries_SP_by_site.png"), p_sites, width = 8, height = 4, dpi = 150)
print(p_sites)

cat("SP sites:", nlevels(sp$site_id), "\n")
cat("Date range:", as.character(min(sp_monthly$date)), "to", as.character(max(sp_monthly$date)), "\n")
cat("Saved: VL_timeseries_SP_monthly.csv, VL_timeseries_SP_total.png, VL_timeseries_SP_by_site.png\n")
