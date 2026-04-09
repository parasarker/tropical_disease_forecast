# Project Data Files

#### ibge_population_byState_2007-2025.csv

-   Variables: "state", "year", "pop_est", "observed"

-   Missing years are filled by linear interpolation (where observed = FALSE).

-   Source: IBGE population estimates

    <https://ftp.ibge.gov.br/Estimativas_de_Populacao/>

#### worldclim_monthly_sites_2007-2024

-   Variables: "site", "year", "month", "tmin_c", "tmax_c", "prec_mm"

-   Source: WorldClim

    <https://geodata.ucdavis.edu/climate/worldclim/2_1/hist/cts4.09/>

#### Met_data temp dir

-   Stage2: forecast met, semi-daily per site (2025-02-10 to 2026-03-13)

    -   Variables: "APCP", "CFRZR", "CICEP", "CRAIN", "CSNOW", "DLWRF", "DSWRF", "ICETK", "LHTFL", "PRES", "PWAT", "RH", "SHTFL", "SNOD", "SOILW", "TCDC", "TMAX", "TMIN", "TMP", "TSOIL", "UGRD", "ULWRF", "USWRF", "VGRD", "WEASD"

-   Stage3: current/analysis met, monthly-ish per site (2026-03-19 and 2026-04-01)

    -   Variables: "air_pressure", "air_temperature", "eastward_wind", "northward_wind", "precipitation_flux", "relative_humidity", "surface_downwelling_longwave_flux_in_air", "surface_downwelling_shortwave_flux_in_air"

-   Source: Met data copied from the SCC until we have an S3 bucket to point to

    [/projectnb/dietzelab/cwebb16/FRP/drivers/stage2_temp](https://scc-ondemand1.bu.edu/pun/sys/dashboard/files/fs//projectnb/dietzelab/cwebb16/FRP/drivers/stage2_temp)

    [/projectnb/dietzelab/cwebb16/FRP/drivers/stage3_temp](https://scc-ondemand1.bu.edu/pun/sys/dashboard/files/fs//projectnb/dietzelab/cwebb16/FRP/drivers/stage2_temp){.uri}
