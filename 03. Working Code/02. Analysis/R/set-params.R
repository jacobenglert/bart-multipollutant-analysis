# Program Name: set-params.R
# Author: Jacob Englert
# Date: 10 January 2024
# Description: Set parameters for BART multipollutant analysis models

# Load Packages -----------------------------------------------------------
library(tidyverse)


# Specify Parameters ------------------------------------------------------

# Run Parameters
outcome <- c("'RESP'","'RESP1'","'ASTHMA'","'ASTHMA1'")

# Model Parameters
seed <- 1
m <- 200
k <- 2
base <- 0.95
power <- 2
num_iter <- 5000
num_burn <- 2500
num_thin <- 5
c <- 0.1
d <- 0.1

# Compile
params <- crossing(outcome, 
                   seed, m, k, base, power, 
                   num_iter, num_burn, num_thin,
                   c, d) |>
  mutate(ID = row_number()) |>
  select(ID, everything())

write_csv(params, here::here('03. Working Code','02. Analysis','Params','params.csv'))
