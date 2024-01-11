# Program Name: fit-model
# Author: Jacob Englert
# Date: 13 November 2023


# Load Packages -----------------------------------------------------------
library(tidyverse)

# Load Data ---------------------------------------------------------------
data <- read_rds(paste0(here::here('02. Analytic Data Set Creation','02. Final Analytic Data Set','atl_data_05-07.rds')))

# Prepare Data ------------------------------------------------------------

# Exclude ZIP codes with no population and those which are not in exposure territory
data <- data[!is.na(data$Tmax),]
data <- data[which(data$POPULATION > 0),]

# Outcome
y <- data$RESP

# Fixed predictors
time_spline <- model.matrix(~splines::ns(yday(data$DATE), df = 7):factor(data$YEAR) - 1)
x1 <- cbind(data$HOLIDAY_FO, time_spline)

# BART predictors
x2 <- select(data, PM25, NO2, O3, CO, Tmax) |> as.matrix()
# write_rds(cor(x2), here::here('03. Working Code','02. Analysis','obs_cor.rds'))

# Geography
geo <- data$geometry

# Offset
offset <- log(data$POPULATION)

# Space-time indicators
s <- as.numeric(as.factor(data$ZIP))
t <- as.numeric(as.factor(data$DATE))

# Remove original data object
rm(data)

# Fit model ---------------------------------------------------------------
source('https://raw.githubusercontent.com/jacobenglert/nb-bart/main/spanbbart.R?token=GHSAT0AAAAAACGFOBGHX2XPYAUAGYCLZ4FYZL4ZPGA')
fit <- spanbbart(x1 = x1, x2 = x2, y = y, s = s, t = t, geo = geo, offset = offset,
                 m = 200, light_store = TRUE,
                 num_iter = 5000, num_burn = 2500, num_thin = 5, seed = 1)
# write_rds(fit, here::here('03. Working Code','02. Analysis','asthma1_05-07.rds'))
fit <- read_rds(paste0(dir, '03. Working Code/02. Analysis/resp_05-07.rds'))


# Examine model fit -------------------------------------------------------

# Variable inclusion proportions
colMeans(fit$var_counts / rowSums(fit$var_counts), na.rm = TRUE)

# Scalar parameters
fit[c('alpha','xi','rho','tau2')] |>
  bind_rows() |>
  mutate(iter = row_number()) |>
  pivot_longer(cols = -iter, names_to = 'parm', values_to = 'value') |>
  mutate(pt = mean(value), 
         l95 = quantile(value, 0.025),
         u95 = quantile(value, 0.975),
         .by = parm) |>
  ggplot(aes(x = iter, y = value)) +
  geom_line() +
  geom_hline(aes(yintercept = pt), color = 'red', lty = 2) +
  geom_hline(aes(yintercept = l95), color = 'red', lty = 2) +
  geom_hline(aes(yintercept = u95), color = 'red', lty = 2) +
  facet_wrap(~parm, scales = 'free', ncol = 1) +
  theme_bw()

plot(log(fit$xi) + fit$alpha, type = 'l', xlab = 'iter', ylab = 'log-mean')

# Random selection of spatial random effects
fit$nu |>
  data.frame() |>
  mutate(iter = row_number()) |>
  pivot_longer(cols = -iter, names_to = 'parm', values_to = 'value') |>
  filter(parm %in% sample(unique(parm), 16)) |>
  mutate(pt = mean(value), 
         l95 = quantile(value, 0.025),
         u95 = quantile(value, 0.975),
         .by = parm) |>
  ggplot(aes(x = iter, y = value)) +
  geom_line() +
  geom_hline(aes(yintercept = pt), color = 'red', lty = 2) +
  geom_hline(aes(yintercept = l95), color = 'red', lty = 2) +
  geom_hline(aes(yintercept = u95), color = 'red', lty = 2) +
  facet_wrap(~parm, scales = 'free') +
  theme_bw()

# Fixed effect coefficients
fit$beta |>
  data.frame() |>
  mutate(iter = row_number()) |>
  pivot_longer(cols = -iter, names_to = 'parm', values_to = 'value') |>
  mutate(pt = mean(value), 
         l95 = quantile(value, 0.025),
         u95 = quantile(value, 0.975),
         .by = parm) |>
  ggplot(aes(x = iter, y = value)) +
  geom_line() +
  geom_hline(aes(yintercept = pt), color = 'red', lty = 2) +
  geom_hline(aes(yintercept = l95), color = 'red', lty = 2) +
  geom_hline(aes(yintercept = u95), color = 'red', lty = 2) +
  facet_wrap(~parm, scales = 'free') +
  theme_bw()



# Accumulated Local Effects -----------------------------------------------

# Custom BART prediction function (returns entire posterior distribution)
pred_fun <- function(model, newdata) model$predict(newdata)

# Load custom Bayesian ALE function
source('https://raw.githubusercontent.com/jacobenglert/nb-bart/main/mcmc_ale.R?token=GHSAT0AAAAAACGFOBGHD7PDCTLKCX6QQS7IZL43KAA')

# Compute first-order ALE for each predictor
K <- 40
ale1 <- parallel::mclapply(1:ncol(x2), \(j) mcmc_ale(x2, fit$bart, pred_fun, j, K), mc.cores = 5)
names(ale1) <- colnames(x2)

# Create data frame of first-order ALEs
ale1df <- mapply(\(x, n, f){
  df <- data.frame(var = n, x = x$x, est = x$est, lcl = x$lcl, ucl = x$ucl)
  return(df)
}, x = ale1, n = names(ale1), SIMPLIFY = FALSE) |>
  bind_rows()

# ale1df <- ale1df |>
#   filter(x != min(x) & x != max(x), .by = var)

# Plot first-order ALEs
ale1df |>
  ggplot(aes(x = x, ymin = exp(lcl), ymax = exp(ucl))) +
  geom_line(aes(y = exp(est))) +
  geom_ribbon(alpha = 0.4) +
  geom_hline(yintercept = 1, lty = 2) +
  facet_wrap(~factor(var, levels = colnames(x2)), scales = 'free') +
  theme_bw() +
  labs(title = 'ALE Plots of First-Order Main Effects',
       x = 'Predictor Value',
       y = 'ALE')


# Compute second-order ALEs for first 5 predictors
pairs <- combn(1:ncol(x2), 2, simplify = FALSE)
ale2 <- parallel::mclapply(pairs, \(j) mcmc_ale(x2, fit$bart, pred_fun, j, K), mc.cores = 2)


# Create data frame of second-order ALEs
ale2df <- mapply(\(x, idx){
  
  n <- colnames(x2)[idx]
  
  x1 <- x$x[[1]]
  x2 <- x$x[[2]]
  
  w1 <- c(diff(x1)[1], diff(x1)) / 2
  w2 <- c(rev(abs(diff(rev(x1)))), rev(abs(diff(x1)))[1]) / 2
  h1 <- c(diff(x2)[1], diff(x2)) / 2
  h2 <- c(rev(abs(diff(rev(x2)))), rev(abs(diff(x2)))[1]) / 2
  
  # Estimate
  f <- as.numeric(x$est)
  
  f2 <- as.numeric(sweep(sweep(x$est, 1, ale1[[idx[1]]]$est, '+'), 2, ale1[[idx[2]]]$est, '+'))
  
  df <- expand.grid(list(x1 = x1, x2 = x2)) |>
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
ale2df$var1 <- factor(ale2df$var1, levels = colnames(x2))
ale2df$var2 <- factor(ale2df$var2, levels = colnames(x2))
  
# ale2df <- ale2df |>
#   filter(x1 != min(x1) & x1 != max(x1) & x2 != min(x2) & x2 != max(x2), .by = c(var1, var2)) |>
#   filter(x1 != min(x1) & x1 != max(x1) & x2 != min(x2) & x2 != max(x2), .by = c(var1, var2))


# Estimate (without Main Effects)
ale2df |>
  ggplot(aes(xmin = x1 - w1, xmax = x1 + w2, ymin = x2 - h1, ymax = x2 + h2)) +
  geom_rect(aes(fill = est)) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
  facet_grid(var2 ~ var1, scales = 'free') +
  theme_bw() +
  labs(title = 'Estimated Second-Order ALEs (excluding Main Effects)',
       x = 'Predictor 1 Value',
       y = 'Predictor 2 Value',
       fill = 'ALE')

# Estimate (with Main Effects)
ale2df |>
  ggplot(aes(xmin = x1 - w1, xmax = x1 + w2, ymin = x2 - h1, ymax = x2 + h2)) +
  geom_rect(aes(fill = est_main)) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
  facet_grid(var2 ~ var1, scales = 'free') +
  theme_bw() +
  labs(title = 'Estimated Second-Order ALEs (including Main Effects)',
       x = 'Predictor 1 Value',
       y = 'Predictor 2 Value',
       fill = 'ALE')







# Partial Dependence ------------------------------------------------------

# # Marginal
# par(mfrow = c(2, 5))
# start <- Sys.time()
# pd <- dbarts::pdbart(fit$bart, xind = colnames(x2), pl = FALSE)
# print(Sys.time() - start)
# plot(pd)
# 
# 
# 
# pd1d <- function(fit, data, var, vals, mc.cores = 1){
#   l_vals <- list(vals)
#   names(l_vals) <- var
#   d1 <- merge(l_vals, data[colnames(data) %in% setdiff(colnames(fit$var_counts), var)])
#   
#   # split_id <- rep(1:nrow(data), each = nrow(d1) / nrow(data))
#   # 
#   # preds <- split(d1, split_id) |>
#   #   parallel::mclapply(\(x) rowMeans(fit$bart$predict(as.matrix(x))), mc.cores = mc.cores) |>
#   #   unlist(use.names = FALSE)
#   # 
#   # # Compute partial dependence
#   # pd <- data.frame(d1, f = preds - mean(fit$alpha)) |>
#   #   summarise(pd = mean(f), .by = c(all_of(var)))
#   
#   
#   preds <- split(d1, d1[[var]]) |>
#     parallel::mclapply(\(x) colMeans(fit$bart$predict(as.matrix(x))), mc.cores = mc.cores) |>
#     unlist(use.names = FALSE)
#   
#   pd <- data.frame(rep(l_vals[[var]], each = length(fit$alpha)), preds - mean(fit$alpha))
#   colnames(pd) <- c(var, 'f')
#   
#   pd <- pd |> 
#     summarise(pd = mean(f),
#               lower = quantile(f, 0.025),
#               upper = quantile(f, 0.975),
#               .by = c(all_of(var)))
#   
#   return(pd)
# }
# 
# var <- 'Tmax'
# vals <- quantile(data$Tmax, probs = seq(0, 1, length.out = 10)) #seq(min(data$Tmax), max(data$Tmax), length.out = 10)
# pd_Tmax <- pd1d(fit, data, var, vals, mc.cores = 3)
# 
# plot(pd_Tmax$pd ~ pd_Tmax$Tmax)
# pd_Tmax |>
#   ggplot(aes(x = Tmax, y = pd, ymin = lower, ymax = upper)) +
#   geom_pointrange()



# Old partial dependence
# Bivariate
par(mfrow = c(1, 1))
start <- Sys.time()
pd2 <- dbarts::pd2bart(fit$bart, xind = O3 + Tmax, pl = FALSE)
print(Sys.time() - start)
plot(pd2)

# Define prediction grid
tmax_vals <- seq(min(data$Tmax), max(data$Tmax), length.out = 10)
pm25_vals <- seq(min(data$PM25), max(data$PM25), length.out = 10)

x2_pdp  <- merge(expand.grid(list(tmax_vals, pm25_vals)), x2[,!colnames(x2) %in% c('Tmax','PM25')])
colnames(x2_pdp) <- c('Tmax','PM25', colnames(x2)[!colnames(x2) %in% c('Tmax','PM25')])

# Make predictions (add back average spatial random effect)
# start <- Sys.time()
# preds <- rowMeans(fit$bart$predict(x2_pdp))
# print(Sys.time() - start)

library(parallel)
start <- Sys.time()
preds <- x2_pdp |> split(rep(1:nrow(x2), each = nrow(x2_pdp) / nrow(x2))) |>
  mclapply(\(x) rowMeans(fit$bart$predict(as.matrix(x))), mc.cores = 4) |>
  unlist(use.names = FALSE)
print(Sys.time() - start)

# Compute partial dependence
pd <- data.frame(x2_pdp, f = preds) |>
  summarise(pd = mean(f), .by = c(Tmax, PM25))

# Estimated bivariate partial dependence surface
bipd_est <- pd |> 
  ggplot(aes(x = Tmax, y = PM25, fill = pd)) +
  geom_tile(col = 'lightgray') +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = mean(fit$alpha)) +
  coord_fixed() +
  theme_bw() + 
  geom_raster(interpolate = TRUE) +
  labs(title = 'Estimated Risk Surface')
bipd_est
