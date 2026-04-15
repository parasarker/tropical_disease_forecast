library(ggplot2)
library(dplyr)
library(tidyr)

# Load disease, met, and population for target sites
monthly_data <- function(
  site_ids = c(150210L, 500270L, 230440L, 170210L, 261110L, 355030L, 312770L, 310620L),
  start_date = as.Date("2007-01-01"), 
  end_date = as.Date("2024-08-01"),
  data_dir = "data",
  disease_url = "https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-read/challenges/project_id=bu4cast/targets/tropical-disease-targets.csv"
) {
  # Disease: monthly case counts
  disease_monthly <- read.csv(disease_url) %>%
    mutate(month = as.Date(substr(datetime, 1, 10))) %>%
    filter(
      site_id %in% site_ids,
      month >= start_date,
      month <= end_date
    ) %>%
    select(site_id, month, observation)

  # Met: WorldClim monthly tmin/tmax/precip
  met_monthly <- read.csv(file.path(data_dir, "worldclim_monthly_sites_2007_2024.csv")) %>%
    mutate(month = as.Date(date)) %>%
    filter(
      site_id %in% site_ids,
      month >= start_date,
      month <= end_date
    ) %>%
    select(site_id, month, tmin_c, tmax_c, prec_mm)

  # Pop: IBGE municipality code (CD_MUN) per forecast site_id (centroids ID / 10).
  cent <- read.csv(file.path(data_dir, "centroids.csv")) %>%
    mutate(site_id = ID %/% 10L)

  site_map <- cent %>%
    transmute(site_id = site_id, CD_MUN = ID) %>%
    distinct() %>%
    filter(site_id %in% site_ids)

  # IBGE municipal population: no 2023
  pop_muni <- read.csv(file.path(data_dir, "municipality_pop_2007_2025.csv")) %>%
    select(CD_MUN, year, pop_est)

  y0 <- as.integer(format(start_date, "%Y"))
  y1 <- as.integer(format(end_date, "%Y"))

  # Annual population time series per site
  pop_site <- site_map %>%
    left_join(pop_muni, by = "CD_MUN") %>%
    filter(year >= y0, year <= y1) %>%
    mutate(year_date = as.Date(paste0(year, "-01-01")))

  # One CD_MUN row per site
  cent_site <- cent %>%
    filter(site_id %in% site_ids) %>%
    distinct(site_id, ID) %>%
    transmute(site_id = site_id, CD_MUN = ID)

  # Site-month table: disease + met + pop (CD_MUN + year)
  monthly_merged <- disease_monthly %>%
    inner_join(met_monthly, by = c("site_id", "month")) %>%
    mutate(year = as.integer(format(month, "%Y"))) %>%
    left_join(cent_site, by = "site_id") %>%
    left_join(pop_muni, by = c("CD_MUN", "year")) %>%
    mutate(pop_scaled = pop_est / 1e6) %>%
    arrange(site_id, month)

  list(
    site_ids = site_ids,
    start_date = start_date,
    end_date = end_date,
    disease_monthly = disease_monthly,
    met_monthly = met_monthly,
    pop_site = pop_site,
    monthly_merged = monthly_merged
  )
}

# Takes in the list from monthly_data()
monthly_plots <- function(d) {
  list(
    disease = ggplot(
      d$disease_monthly,
      aes(x = month, y = observation, color = factor(site_id), group = site_id)
    ) +
      geom_line(alpha = 0.8, linewidth = 0.7) +
      labs(
        title = "Monthly Disease Observations by Site",
        x = "Month",
        y = "Cases",
        color = "Site ID"
      ) +
      theme_minimal(),
    met = ggplot(
      d$met_monthly %>%
        pivot_longer(
          c(tmin_c, tmax_c, prec_mm),
          names_to = "variable",
          values_to = "value"
        ),
      aes(x = month, y = value, color = factor(site_id), group = site_id)
    ) +
      geom_line(alpha = 0.8, linewidth = 0.6) +
      facet_wrap(~variable, scales = "free_y", ncol = 1) +
      labs(
        title = "WorldClim Monthly Met Variables by Site",
        x = "Month",
        y = NULL,
        color = "Site ID"
      ) +
      theme_minimal(),
    pop = ggplot(
      d$pop_site,
      aes(x = year_date, y = pop_est, color = factor(site_id), group = site_id)
    ) +
      geom_line(alpha = 0.85, linewidth = 0.8) +
      geom_point(alpha = 0.7, size = 1.2) +
      labs(
        title = "Municipality Population by Site (Annual)",
        x = "Year",
        y = "Population",
        color = "Site ID"
      ) +
      theme_minimal()
  )
}

# d <- monthly_data()
# p <- monthly_plots(d)
# print(p$disease); print(p$met); print(p$pop)
