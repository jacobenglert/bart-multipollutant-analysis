dir <- "/Volumes/rsphprojects-ts/envision/Analysis/Jacob/BART/"
data <- read_rds(paste0(dir, '02. Analytic Data Set Creation/02. Final Analytic Data Set/atl_data_05-07.rds'))

shp <- data |>
  filter(POPULATION > 0) |>
  filter(DATE == min(DATE), .by = ZIP)
colnames(data)

shp |> ggplot(aes(geometry=geometry)) + geom_sf()

shp |> select(POPULATION, ZIP, geometry) |> write_rds(here::here('Data','Clean','GAZIP05.rds'))
