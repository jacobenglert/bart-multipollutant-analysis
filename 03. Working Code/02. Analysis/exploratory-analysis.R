# Program Name: fit-model
# Author: Jacob Englert
# Date: 13 November 2023


# Load Packages -----------------------------------------------------------
library(tidyverse)


# Load Data ---------------------------------------------------------------
dir <- "/Volumes/rsphprojects-ts/envision/Analysis/Jacob/BART/"
data <- read_rds(paste0(dir, '02. Analytic Data Set Creation/02. Final Analytic Data Set/atl_data_05-07.rds'))
data <- na.omit(data)
data <- filter(data, POPULATION > 0)



# Exploratory Analysis ----------------------------------------------------
library(splines)

# data$TimeSpline <- ns(yday(data$DATE), df = 7) * data$YEAR

# Try fitting negative binomial glm
fit.nb.resp <- MASS::glm.nb(RESP ~ ns(yday(DATE), df = 7) : as.factor(YEAR) + 
                              PM25 + O3 + NO2 + Tmax + offset(log(POPULATION)), 
                            data = data)
summary(fit.nb.resp)$coefficients[2:5,]

fit.nb.resp1 <- MASS::glm.nb(RESP1 ~ ns(yday(DATE), df = 7) : as.factor(YEAR) + 
                               PM25 + O3 + NO2 + Tmax + offset(log(POPULATION)), 
                             data = data)
summary(fit.nb.resp1)$coefficients[2:5,]

fit.nb.asth <- MASS::glm.nb(ASTHMA ~ ns(yday(DATE), df = 7) : as.factor(YEAR) + 
                              PM25 + O3 + NO2 + Tmax + offset(log(POPULATION)), 
                            data = data)
summary(fit.nb.asth)$coefficients[2:5,]

fit.nb.asth1 <- MASS::glm.nb(ASTHMA1 ~ ns(yday(DATE), df = 7) : as.factor(YEAR) + 
                               PM25 + O3 + NO2 + Tmax + offset(log(POPULATION)), 
                             data = data)
summary(fit.nb.asth1)$coefficients[2:5,]


# # Try fitting model with spatial random effects
# library(lme4)
# test <- glmer(RESP ~ (1 | ZIP) + ns(yday(DATE), df = 7) : as.factor(YEAR) + 
#                 PM25 + O3 + NO2 + Tmax + offset(log(POPULATION)), 
#               data = data, family = 'poisson',
#               control = glmerControl(optimizer = 'bobyqa'), nAGQ = 50)
# summary(test)
# 
# # Try fitting zero-inflated model
# library(glmmTMB)
# test <- glmmTMB(RESP ~ (1 | ZIP) + ns(yday(DATE), df = 7) : as.factor(YEAR) + 
#                   PM25 + O3 + NO2 + Tmax + offset(log(POPULATION)), 
#                 data = data, family = nbinom2(), 
#                 ziformula = ~1)