# Program Name: set-analysis-params.R
# Author: Jacob Englert
# Date: 13 February 2024
# Description: Set parameters for BART multipollutant analysis models

# Load Packages -----------------------------------------------------------
library(tidyverse)


# Specify Parameters ------------------------------------------------------

# Outcome of interest
outcome <- c("ASTHMA1")

# Exposures of interest
exposures <- c("PM25","NO2","O3","CO","Tmax")
exposure_combos <- list()
for (i in 1:length(exposures)) {
  exposure_combos <- c(exposure_combos, combn(exposures, i, simplify = FALSE))
}

# Model Parameters
seed <- 1
m <- 50
k <- 2
base <- 0.95
power <- 2
num_iter <- 10000
num_burn <- 5000
num_thin <- 10
c <- 0.1
d <- 0.1
sparse <- TRUE
soft <- TRUE
K <- 40


# Compile
params <- crossing(outcome, exposures = exposure_combos,
                   seed, m, k, base, power, 
                   num_iter, num_burn, num_thin,
                   c, d, sparse, soft) |>
  mutate(key = row_number()) |>
  select(key, everything())


# Export
write_rds(params, here::here('03. Working Code','02. Analysis','Params','analysis-params.rds'))
