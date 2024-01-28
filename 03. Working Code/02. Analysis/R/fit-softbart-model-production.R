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
x2 <- SoftBart::quantile_normalize_bart(x2)

# Space-time indicators
s <- as.numeric(as.factor(data$ZIP))
t <- as.numeric(as.factor(data$DATE))


# Fit Model ---------------------------------------------------------------
m <- 50
fit <- soft_spanbbart(x1 = x1, x2 = x2, y = y, s = s, t = t, 
                 geo = data$geometry, offset = log(data$POPULATION),
                 seed = seed, m = m, k = k, base = base, power = power,
                 c = c, d = d,
                 num_iter = num_iter, num_burn = num_burn, num_thin = num_thin,
                 light_store = TRUE)

# Export Model Fit
write_rds(fit, here::here('03. Working Code','02. Analysis','Results', paste0('softbart-fit-', sprintf('%02d', index), '.rds')))



# Compute ALE -------------------------------------------------------------

# NEED TO DO THIS PART IN THE SAME SESSION AS THE MODEL IS FIT!!!

# Prediction function
pred_fun <- function(model, newdata) {
  sapply(seq.int(num_burn + 1, num_iter, num_thin), \(i) model$predict_iteration(newdata, i))
}


K <- 100
var_names <- colnames(x2)
P <- ncol(fit$var_counts)


# Compute First-Order ALE -------------------------------------------------

ale1 <- parallel::mclapply(1:P, \(j) mcmc_ale(x2, fit$bart, pred_fun, j, K), mc.cores = P)
names(ale1) <- var_names

# Create data frame of first-order ALEs
ale1df <- mapply(\(x, n){
  df <- data.frame(var = n, x = x$x, est = x$est, lcl = x$lcl, ucl = x$ucl)
  return(df)
}, x = ale1, n = names(ale1), SIMPLIFY = FALSE) |>
  bind_rows()
ale1df$var <- factor(ale1df$var, levels = var_names)


# Compute Second-Order ALE ------------------------------------------------

pairs <- combn(1:P, 2, simplify = FALSE)
ale2 <- parallel::mclapply(pairs, \(j) mcmc_ale(x2, fit$bart, pred_fun, j, K), mc.cores = P)

# Create data frame of second-order ALEs
ale2df <- mapply(\(x, idx){
  
  n <- var_names[idx]
  
  x1 <- x$x[[1]]
  x2 <- x$x[[2]]
  
  w1 <- c(diff(x1)[1], diff(x1)) / 2
  w2 <- c(rev(abs(diff(rev(x1)))), rev(abs(diff(x1)))[1]) / 2
  h1 <- c(diff(x2)[1], diff(x2)) / 2
  h2 <- c(rev(abs(diff(rev(x2)))), rev(abs(diff(x2)))[1]) / 2
  
  # Estimate
  f <- as.numeric(x$est)
  f2 <- as.numeric(sweep(sweep(x$est, 1, ale1[[idx[1]]]$est, '+'), 2, ale1[[idx[2]]]$est, '+'))
  
  df <- list(x1 = x1, x2 = x2) |>
    expand.grid() |>
    mutate(w1 = rep(w1, times = length(x$x[[2]])),
           w2 = rep(w2, times = length(x$x[[2]])),
           h1 = rep(h1, each = length(x$x[[1]])),
           h2 = rep(h2, each = length(x$x[[1]]))) |>
    mutate(est = f,
           lcl = as.numeric(x$lcl),
           ucl = as.numeric(x$ucl),
           est_main = f2) |>
    mutate(var1 = n[1], var2 = n[2])
  
  return(df)
}, x = ale2, idx = pairs, SIMPLIFY = FALSE) |>
  bind_rows()
ale2df$var1 <- factor(ale2df$var1, levels = var_names)
ale2df$var2 <- factor(ale2df$var2, levels = var_names)


# Output Results ----------------------------------------------------------
write_rds(list(ale1 = ale1df, ale2 = ale2df), 
          here::here('03. Working Code','02. Analysis','Results', paste0('softbart-fit-', sprintf('%02d', index), '-ale-', K ,'.rds')))


