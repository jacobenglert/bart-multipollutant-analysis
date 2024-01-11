# Program Name: compute-ale
# Author: Jacob Englert
# Date: 11 January 2024
# Purpose: compute first- and second-order ALEs for fitted models

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(nbbart)


# Get Model Fit Index -----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
names(args) <- c('index','K')
index <- as.numeric(args['index'])
K <- as.numeric(args['K'])


# Load Model Fit ----------------------------------------------------------

# Model
fit <- read_rds(here::here('03. Working Code','02. Analysis','Results', paste0('fit-', sprintf('%02d', index), '.rds')))

# Prediction function
pred_fun <- function(model, newdata) model$predict(newdata)


# Load Training Data ------------------------------------------------------
x2 <- read_rds(here::here('02. Analytic Data Set Creation','02. Final Analytic Data Set','atl_data_05-07.rds')) |>
  sf::st_drop_geometry() |>
  na.omit() |>
  filter(POPULATION > 0) |>
  select(all_of(colnames(fit$var_counts))) |>
  as.matrix()

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
          here::here('03. Working Code','02. Analysis','Results', paste0('fit-', sprintf('%02d', index), '-ale-', K ,'.rds')))

