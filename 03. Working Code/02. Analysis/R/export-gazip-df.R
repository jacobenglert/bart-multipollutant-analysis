
# Export ATL Zip code shapefile for use in simulations

library(tidyverse)

data <- read_rds(here::here('02. Analytic Data Set Creation','02. Final Analytic Data Set','atl_data_05-07.rds'))

# Only use populated regions
shp <- data |>
  filter(POPULATION > 0) |>
  filter(DATE == min(DATE), .by = ZIP) # Only keep one record for each ZIP

# Visualize regions
shp |> 
  ggplot(aes(geometry = geometry)) + 
  geom_sf()

# Export
shp |> 
  select(POPULATION, ZIP, geometry) |> 
  write_rds(here::here('03. Working Code','02. Analysis','Misc','GAZIP05.rds'))
