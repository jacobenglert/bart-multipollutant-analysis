# Program Name: run-analysis.R
# Author: Jacob Englert
# Date: 13 February 2024
# Purpose: Fit spanbbart model to ATL data

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(splines)
library(nbbart)
library(sf)


# Get Parameters ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
key <- args[1]
params <- read_rds(here::here('03. Working Code','02. Analysis','Params','analysis-params.rds'))
for(i in 1:ncol(params)){
  name <- names(params)[i]
  value <- params[key, i, drop = TRUE]
  if (is.numeric(value)) assign(name, eval(parse(text = value)))
  else if (is.character(value)) assign(name, value)
  else if (is.logical(value)) assign(name, value)
  else if (is.list(value)) assign(name, unlist(value))
}


# Load Data ---------------------------------------------------------------
data <- read_rds(here::here('02. Analytic Data Set Creation','02. Final Analytic Data Set','atl_data_05-07.rds'))

# Compute 3 day moving average
data <- data |>
  arrange(ZIP, DATE) |>
  mutate(across(all_of(exposures), ~zoo::rollmean(x = ., k = 3, fill = NA, align = 'right')), .by = c(ZIP,YEAR))
data <- na.omit(data) # drop values with missing data (typically exposures)


# Prepare Data ------------------------------------------------------------

# Outcome
y <- data[[outcome]]

# Design matrices
time_spline <- model.matrix(~ns(yday(data$DATE), df = 7):factor(data$YEAR) - 1)
x1 <- cbind(data$HOLIDAY_FO, time_spline)
x2 <- data[exposures] |> as.matrix()
x2_scaled <- apply(x2, 2, \(x) (x - min(x))/ (max(x) - min(x))) # Scale BART predictors for SoftBART model

# Space-time indicators
s <- as.numeric(as.factor(data$ZIP))
t <- as.numeric(as.factor(data$DATE))


# Fit Model ---------------------------------------------------------------
fit <- soft_spanbbart(x1 = x1, x2 = x2_scaled, y = y, s = s, t = t, 
                      geo = data$geometry, offset = log(data$POPULATION),
                      seed = seed, m = m, k = k, base = base, power = power,
                      c = c, d = d, sparse = sparse, soft = soft,
                      num_iter = num_iter, num_burn = num_burn, num_thin = num_thin,
                      light_store = TRUE)


# Compute ALE -------------------------------------------------------------

# NEED TO DO THIS PART IN THE SAME SESSION AS THE MODEL IS FIT!!!

# Prediction function
pred_fun <- function(model, newdata) {
  sapply(seq.int(num_burn + 1, num_iter, num_thin), \(i) model$predict_iteration(newdata, i))
}

var_names <- colnames(x2)
P <- ncol(fit$var_counts)

# Rescale x values back to original scale
rescale1 <- function (x2, ale) {
  xmin <- min(x2[,ale$var])
  xmax <- max(x2[,ale$var])
  ale$x <- (ale$x * (xmax - xmin)) + xmin
  return (ale)
}

# Rescale x values and plotting coordinates back to original scale
rescale2 <- function (x2, ale) {
  
  x1old <- unique(ale$x1)
  x1min <- min(x2[,ale$var1])
  x1max <- max(x2[,ale$var1])
  x1new <- (x1old * (x1max - x1min)) + x1min
  
  x2old <- unique(ale$x2)
  x2min <- min(x2[,ale$var2])
  x2max <- max(x2[,ale$var2])
  x2new <- (x2old * (x2max - x2min)) + x2min
  
  w1new <- c(diff(x1new)[1], diff(x1new)) / 2
  w2new <- c(rev(abs(diff(rev(x1new)))), rev(abs(diff(x1new)))[1]) / 2
  h1new <- c(diff(x2new)[1], diff(x2new)) / 2
  h2new <- c(rev(abs(diff(rev(x2new)))), rev(abs(diff(x2new)))[1]) / 2

  index <- 1
  for (i in 1:length(x1new)) {
    for (j in 1:length(x2new)) {
      ale$x1[index] <- x1new[i]
      ale$x2[index] <- x2new[j]
      ale$w1[index] <- w1new[i]
      ale$w2[index] <- w2new[i]
      ale$h1[index] <- h1new[j]
      ale$h2[index] <- h2new[j]
      index <- index + 1
    }
  }
  return (ale)
}

# Compute First-Order ALE -------------------------------------------------

ale1 <- parallel::mclapply(1:P, \(j) mcmc_ale(x2_scaled, fit$bart, pred_fun, j, K, center = 'avg_pred'), mc.cores = P)

# Create data frame of first-order ALEs
ale1df <- lapply(ale1, \(ale) rescale1(x2, ale)) |>
  do.call(what = rbind)
ale1df$var <- factor(ale1df$var, levels = var_names)


# Compute Second-Order ALE ------------------------------------------------

pairs <- combn(1:P, 2, simplify = FALSE)
ale2 <- parallel::mclapply(pairs, \(j) mcmc_ale(x2_scaled, fit$bart, pred_fun, 
                                                j, K, center = 'avg_pred'), 
                           mc.cores = P)

ale2df <- lapply(ale2, \(ale) rescale2(x2, ale)) |>
  do.call(what = rbind)
ale2df$var1 <- factor(ale2df$var1, levels = var_names)
ale2df$var2 <- factor(ale2df$var2, levels = var_names)


# Compute Second-Order ALE with Main Effects ------------------------------

ale3 <- parallel::mclapply(pairs, \(j) mcmc_ale(x2_scaled, fit$bart, pred_fun, 
                                                j, K, center = 'avg_pred',
                                                include_main_effects = TRUE),
                           mc.cores = P)

ale3df <- lapply(ale3, \(ale) rescale2(x2, ale)) |>
  do.call(what = rbind)
ale3df$var1 <- factor(ale3df$var1, levels = var_names)
ale3df$var2 <- factor(ale3df$var2, levels = var_names)


# Output Results ----------------------------------------------------------
write_rds(list(key = key, fit = fit, ale1 = ale1df, ale2 = ale2df, ale3 = ale3df), 
          here::here('03. Working Code','02. Analysis','Results', paste0(key, '.rds')))


