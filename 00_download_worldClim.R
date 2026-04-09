# Pull WorldClim monthly tmin/tmax/prec for our 8 sites (2007-2024)

library(terra)

# selected sites
site_ids<- c(150210, 500270, 230440, 170210, 261110, 355030, 312770, 310620)

res <- "10m" # 10m =coarsest, there's also 5m or 2.5m

base <- "https://geodata.ucdavis.edu/climate/worldclim/2_1/hist/cts4.09/"
zipdir <- "data/worldclim_zips"
tifdir <- "data/worldclim_tifs"
outfile <- "tables/worldclim_monthly_sites_2007_2024.csv"

dir.create(zipdir, FALSE, TRUE)
dir.create(tifdir, FALSE, TRUE)
dir.create(dirname(outfile), FALSE, TRUE)

cent <- read.csv("centroids.csv", stringsAsFactors = FALSE)
cent$site_id <- cent$ID %/% 10L  # site_ids don't have the last digit in centroid_ids
sites <- cent[cent$site_id %in% site_ids, ]
sites <- sites[match(site_ids, sites$site_id), ]
pts <- vect(sites, geom = c("X", "Y"), crs = "EPSG:4326")

# worldclim stores 10 years in each zip file (except last chunk which ends 2024)
which_zip <- function(y) {
  if (y <= 2009) {
    return("2000-2009")
  }
  if (y <= 2019) {
    return("2010-2019")
  }
  "2020-2024"
}

vars <- c("tmin", "tmax", "prec")
ranges <- c("2000-2009", "2010-2019", "2020-2024")

# 3 time chunks per variable (2000-2009, 2010-2019, 2020-2024) -> 9 zip files total
# we don't want to re-download/unzip those same files each rerun
for (vr in vars) {
  for (rg in ranges) {
    zipname <- paste0("wc2.1_cruts4.09_", res, "_", vr, "_", rg, ".zip")
    zippath <- file.path(zipdir, zipname)
    if (!file.exists(zippath)) {
      download.file(paste0(base, zipname), zippath, mode = "wb")
    }
    folder <- file.path(tifdir, tools::file_path_sans_ext(zipname))
    if (!dir.exists(folder) || length(list.files(folder, pattern = "\\.tif$")) == 0) {
      dir.create(folder, FALSE, TRUE)
      unzip(zippath, exdir = folder)
    }
  }
}

# loop every month, grab pixel at each site
months <- seq(as.Date("2007-01-01"), as.Date("2024-12-31"), by = "month")
chunks <- vector("list", length(months))

for (i in seq_along(months)) {
  d <- months[i]
  ym <- format(d, "%Y-%m")
  y <- as.integer(substr(ym, 1, 4))  #chars 1:4 (year) from "YYYY-MM"
  mo <- as.integer(substr(ym, 6, 7)) # chars 6:7 (month)
  rg <- which_zip(y)
  tmin_f <- file.path(
    tifdir,
    paste0("wc2.1_cruts4.09_", res, "_tmin_", rg),
    paste0("wc2.1_cruts4.09_", res, "_tmin_", ym, ".tif")
  )
  tmax_f <- file.path(
    tifdir,
    paste0("wc2.1_cruts4.09_", res, "_tmax_", rg),
    paste0("wc2.1_cruts4.09_", res, "_tmax_", ym, ".tif")
  )
  prec_f <- file.path(
    tifdir,
    paste0("wc2.1_cruts4.09_", res, "_prec_", rg),
    paste0("wc2.1_cruts4.09_", res, "_prec_", ym, ".tif")
  )
  
  ex <- extract(rast(tmin_f), pts, ID = FALSE)
  tmin <- as.numeric(ex[[length(ex)]])
  ex <- extract(rast(tmax_f), pts, ID = FALSE)
  tmax <- as.numeric(ex[[length(ex)]])
  ex <- extract(rast(prec_f), pts, ID = FALSE)
  prec <- as.numeric(ex[[length(ex)]])
  
  # one table per month
  chunks[[i]] <- data.frame(
    site_id = as.integer(sites$site_id),
    name = as.character(sites$name),
    state = as.character(sites$state),
    climreg = as.character(sites$climreg),
    lon = as.numeric(sites$X),
    lat = as.numeric(sites$Y),
    year = y,
    month = mo,
    date = paste0(ym, "-01"),
    tmin_c = tmin,
    tmax_c = tmax,
    prec_mm = prec,
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, chunks)
rownames(out) <- NULL #to get default numbering

# output is a long format table with one row per site-month
write.csv(out, outfile, row.names = FALSE)