# Program Name: create-analysis-dataset
# Author: Jacob Englert
# Date: 13 November 2023


# Load Required Packages --------------------------------------------------
library(tidyverse)
library(tigris)
library(sf)


# Obtain list of Atlanta 5-County area zip codes --------------------------

# Load zip code shapefile
load(here::here('01. Data','Exposure Data','Shapefile','USZIP5_2019_GA.RData'))
zip_codes <- state.shp |>
  mutate(POPULATION = as.numeric(POPULATION)) |>
  filter(POPULATION > 0)
setdiff(state.shp$ZIP, zip_codes$ZIP)
rm(state.shp)

# Load county shapefile
counties <- tigris::counties(state = 'GA', class = "sf", year = 2019)
counties <- counties |>
  filter(NAME %in% c('Fulton','DeKalb','Gwinnett','Cobb','Clayton')) |>
  sf::st_transform('WGS84')

# Identify zip codes in the 5-county ATL
atl_zips_idx <- sf::st_intersects(zip_codes, counties, sparse = FALSE) |>
  apply(1, any)
atl_zips <- zip_codes$ZIP[atl_zips_idx]

# plot(sf::st_geometry(counties), col = 'red')
# plot(sf::st_geometry(zip_codes[which(zip_codes$ZIP %in% atl_zips),]), add = TRUE)



# Load Health Data --------------------------------------------------------
h_data <- read_csv(here::here('01. Data','Health Data','FINAL_ED_ZIP_DAILY_COUNTS.csv'))
h_data_atl <- h_data |>
  filter(PATZIP %in% atl_zips) |>
  mutate(DATE = as.Date(DATE, format = "%m/%d/%Y")) |>
  arrange(PATZIP, DATE)

# Identify duplicates
dupes <- h_data_atl |> mutate(n = max(row_number()), .by = c(DATE, PATZIP)) |> filter(n > 1)

# Add duplicate rows together
h_data_atl <- h_data_atl |>
  group_by(PATZIP, DATE) |>
  summarise(across(RESP:ASTHMA1, sum)) |>
  ungroup()

# Load Exposure Data ------------------------------------------------------

# Air pollution data
e_data1 <- data.frame()
for(i in 2005:2019){
  load(here::here('01. Data','Exposure Data','Air Pollution','georgia', paste0(i, '_georgia.RData')))
  e_data1 <- bind_rows(e_data1, dat.zip)
}
rm(dat.zip)
e_data1_atl <- filter(e_data1, ZIP %in% atl_zips)
length(unique(e_data1_atl$ZIP))
table(table(e_data1_atl$ZIP))
table(e_data1_atl$ZIP)[which(table(e_data1_atl$ZIP) < 5478)]

# Max temperature data
e_data2_atl <- data.frame()
for(i in 2005:2018){
  e_zips <- read_csv(here::here('01. Data','Exposure Data','Met',paste0('ZIP_order_', i, '.zip')),
                     show_col_types = FALSE, progress = FALSE,
                     col_names = FALSE)[[1]]
  atl_e_zips_idx <- which(e_zips %in% atl_zips)
  
  all_dates_in_year <- seq(as.Date(paste0(i, "/01/01")),
                           as.Date(paste0(i, "/12/31")),
                           by = "days")
  
  tmax_data <- read_csv(here::here('01. Data','Exposure Data','Met',paste0('tmax_', i, '_ZIP.zip')), 
                        show_col_types = FALSE, progress = FALSE,
                        col_names = FALSE)[atl_e_zips_idx,] |>
    mutate(ZIP = as.character(e_zips[atl_e_zips_idx])) |>
    pivot_longer(cols = -ZIP, values_to = 'Tmax') |>
    mutate(DATE = rep(all_dates_in_year, times = length(atl_e_zips_idx))) |>
    select(-name)
  
  tmin_data <- read_csv(here::here('01. Data','Exposure Data','Met',paste0('tmin_', i, '_ZIP.zip')), 
                        show_col_types = FALSE, progress = FALSE,
                        col_names = FALSE)[atl_e_zips_idx,] |>
    mutate(ZIP = as.character(e_zips[atl_e_zips_idx])) |>
    pivot_longer(cols = -ZIP, values_to = 'Tmin') |>
    mutate(DATE = rep(all_dates_in_year, times = length(atl_e_zips_idx))) |>
    select(-name)
  
  e_data2_atl <- bind_rows(e_data2_atl, left_join(tmax_data, tmin_data, by = c('ZIP','DATE')))
  rm(tmax_data, tmin_data, all_dates_in_year, atl_e_zips_idx, e_zips)
  print(paste0('Year ', i, ' completed.'))
}

length(unique(e_data2_atl$ZIP))
table(table(e_data2_atl$ZIP))
table(e_data2_atl$ZIP)[which(table(e_data2_atl$ZIP) < 5113)]


# Time Data ---------------------------------------------------------------

t_data <- read_csv(here::here('01. Data','Time Data','timevars_93_20.csv')) |>
  mutate(DATE = as.Date(DATE, format = "%m/%d/%Y")) |>
  select(DATE, MONTH, YEAR, HOLIDAY_FO, WARM_ATL, COLD_ATL)


# Join Data ---------------------------------------------------------------

atl_data <- h_data_atl |>
  left_join(e_data1_atl, by = c('PATZIP' = 'ZIP', 'DATE' = 'date')) |>
  left_join(e_data2_atl, by = c('PATZIP' = 'ZIP', 'DATE')) |>
  left_join(t_data, by = c('DATE')) |>
  left_join(zip_codes, by = c('PATZIP' = 'ZIP')) |>
  select(ZIP = PATZIP, DATE, RESP:ASTHMA1, PM25:Tmin, MONTH:COLD_ATL, POPULATION, geometry) |>
  filter(YEAR %in% 2005:2018)


# Export Data -------------------------------------------------------------

# For now, focus on warm season for 2005-2007
atl_data |>
  filter(WARM_ATL == 1 & YEAR %in% 2005:2007) |>
  write_rds(here::here('02. Analytic Data Set Creation','02. Final Analytic Data Set','atl_data_05-07.rds'))

# All data
atl_data |>
  write_rds(here::here('02. Analytic Data Set Creation','02. Final Analytic Data Set','atl_data.rds'))

# All data without geometries
atl_data |>
  sf::st_drop_geometry() |>
  write_rds(here::here('02. Analytic Data Set Creation','02. Final Analytic Data Set','atl_data_df.rds'))


