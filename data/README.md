# Project Data Files

#### state_pop_by_year.csv

-   Brazil annual population data by state ("state", "year", "pop_est", "observed").

-   Values are aggregated to state level.

-   Missing years are filled by linear interpolation (where observed = FALSE).

-   Source: IBGE population estimates

    <https://ftp.ibge.gov.br/Estimativas_de_Populacao/>

#### Met_data temp dirs

-   Stage 3 and Stage2 met data copied from the SCC until we have an S3 bucket to point to

-   Semi-daily met data by site in Sao Paolo

-   Source:
-     SCC "/projectnb/dietzelab/cwebb16/FRP/drivers/stage3_temp"
-     SCC "/projectnb/dietzelab/cwebb16/FRP/drivers/stage2_temp"
