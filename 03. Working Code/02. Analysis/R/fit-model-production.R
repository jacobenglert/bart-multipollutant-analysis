# Program Name: fit-model-production
# Author: Jacob Englert
# Date: 10 January 2024
# Purpose: Fit spanbbart model to ATL data

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(splines)
library(nbbart)
library(sf)


# Get Parameters ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
names(args) <- c('index')
index <- as.numeric(args['index'])
params <- read_csv(here::here('03. Working Code','02. Analysis','Params','params.csv'), show_col_types = F)[,-1]
for(i in 1:ncol(params)){
  assign(names(params[index,i]), eval(parse(text = params[index, i, drop = TRUE])))
}


# Load Data ---------------------------------------------------------------
data <- read_rds(here::here('02. Analytic Data Set Creation','02. Final Analytic Data Set','atl_data_05-07.rds'))
data <- na.omit(data)
data <- filter(data, POPULATION > 0)


# Prepare Data ------------------------------------------------------------

# Outcome
y <- data[[outcome]]

# Design matrices
time_spline <- model.matrix(~ns(yday(data$DATE), df = 7):factor(data$YEAR) - 1)
x1 <- cbind(data$HOLIDAY_FO, time_spline)
x2 <- dplyr::select(data, PM25, NO2, O3, CO, Tmax) |> as.matrix()

# Space-time indicators
s <- as.numeric(as.factor(data$ZIP))
t <- as.numeric(as.factor(data$DATE))


# Fit Model ---------------------------------------------------------------
fit <- spanbbart(x1 = x1, x2 = x2, y = y, s = s, t = t, 
                 geo = data$geometry, offset = log(data$POPULATION),
                 seed = seed, m = m, k = k, base = base, power = power,
                 c = c, d = d,
                 num_iter = num_iter, num_burn = num_burn, num_thin = num_thin,
                 light_store = TRUE)

write_rds(fit, here::here('03. Working Code','02. Analysis','Results', paste0('fit-', sprintf('%02d', index), '.rds')))


