# Program Name: fit-model
# Author: Jacob Englert
# Date: 13 November 2023


# Load Packages -----------------------------------------------------------
library(tidyverse)
library(splines)
source('spanbbart.R')

# Load Data ---------------------------------------------------------------
dir <- "/Volumes/rsphprojects-ts/envision/Analysis/Jacob/BART/"
data <- read_rds(paste0(dir, '02. Analytic Data Set Creation/02. Final Analytic Data Set/atl_data_05-07.rds'))
data <- na.omit(data)
data <- filter(data, POPULATION > 0)


# Prepare Data ------------------------------------------------------------

# Outcome
y <- data$ASTHMA1
time_spline <- model.matrix(~ns(yday(data$DATE), df = 7):factor(data$YEAR) - 1)

# Design matrices
x1 <- cbind(data$HOLIDAY_FO, time_spline)
x2 <- dplyr::select(data, PM25, NO2, O3, Tmax) |> as.matrix()

# Space-time indicators
s <- as.numeric(as.factor(data$ZIP))
t <- as.numeric(as.factor(data$DATE))


# Fit Model ---------------------------------------------------------------
fit <- spanbbart(x1 = x1, x2 = x2, y = y, s = s, t = t, 
                 geo = data$geometry, offset = log(data$POPULATION),
                 m = 100, light_store = TRUE,
                 num_iter = 5000, num_burn = 2500, num_thin = 5, seed = 1)
sa


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



# ALE Plots ---------------------------------------------------------------
mcmc_pred_func <- function(model, newdata) model$predict(newdata)
ALE.1 <- mcmc_ale(x2, fit$bart, pred_fun = mcmc_pred_func, J = 1, K = 20)
ALE.2 <- mcmc_ale(x2, fit$bart, pred_fun = mcmc_pred_func, J = 2, K = 20)
ALE.3 <- mcmc_ale(x2, fit$bart, pred_fun = mcmc_pred_func, J = 3, K = 20)
ALE.4 <- mcmc_ale(x2, fit$bart, pred_fun = mcmc_pred_func, J = 4, K = 20)

ALE.4.1 <- mcmc_ale(x2, fit$bart, pred_fun = mcmc_pred_func, J = c(4,1), K = 20)
ALE.4.2 <- mcmc_ale(x2, fit$bart, pred_fun = mcmc_pred_func, J = c(4,2), K = 20)
ALE.4.3 <- mcmc_ale(x2, fit$bart, pred_fun = mcmc_pred_func, J = c(4,3), K = 20)


data.frame(x = ALE.1$x.values, y = ALE.1$f.values, l = ALE.1$f.lowers, u = ALE.1$f.uppers) |>
  ggplot(aes(x = x, y = y, ymin = l, ymax = u)) +
  geom_line() +
  geom_ribbon(alpha = 0.4) +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(-1,1))

data.frame(x = ALE.2$x.values, y = ALE.2$f.values, l = ALE.2$f.lowers, u = ALE.2$f.uppers) |>
  ggplot(aes(x = x, y = y, ymin = l, ymax = u)) +
  geom_line() +
  geom_ribbon(alpha = 0.4)

data.frame(x = ALE.3$x.values, y = ALE.3$f.values, l = ALE.3$f.lowers, u = ALE.3$f.uppers) |>
  ggplot(aes(x = x, y = y, ymin = l, ymax = u)) +
  geom_line() +
  geom_ribbon(alpha = 0.4)

data.frame(x = ALE.4$x.values, y = ALE.4$f.values, l = ALE.4$f.lowers, u = ALE.4$f.uppers) |>
  ggplot(aes(x = x, y = y, ymin = l, ymax = u)) +
  geom_line() +
  geom_ribbon(alpha = 0.4)


par(mfrow = c(1,3))
f.41 <- sweep(sweep(ALE.4.1$f.values, 1, ALE.4$f.values, "+"), 2, ALE.1$f.values, "+")
fields::image.plot(ALE.4.1$x.values[[1]], ALE.4.1$x.values[[2]], f.41, 
                   xlab = "Tmax", ylab = "PM2.5", col = hcl.colors(10, palette = 'Blues'))
#contour(ALE.41$x.values[[1]], ALE.41$x.values[[2]], f.41, add = TRUE, drawlabels=TRUE)

f.42 <- sweep(sweep(ALE.4.2$f.values, 1, ALE.4$f.values, "+"), 2, ALE.2$f.values, "+")
fields::image.plot(ALE.4.2$x.values[[1]], ALE.4.2$x.values[[2]], f.42, 
                   xlab = "Tmax", ylab = "NO2", col = hcl.colors(10, palette = 'Blues'))
#contour(ALE.42$x.values[[1]], ALE.42$x.values[[2]], f.42, add = TRUE,drawlabels=TRUE)

f.43 <- sweep(sweep(ALE.4.3$f.values, 1, ALE.4$f.values, "+"), 2, ALE.3$f.values, "+")
fields::image.plot(ALE.4.3$x.values[[1]], ALE.4.3$x.values[[2]], f.43, 
                   xlab = "Tmax", ylab = "O3", col = hcl.colors(10, palette = 'Blues'))
#contour(ALE.43$x.values[[1]], ALE.43$x.values[[2]], f.43, add = TRUE, drawlabels=TRUE)
















library(ALEPlot)
my_pred_func <- function(X.model, newdata) as.numeric(rowMeans(X.model$predict(newdata)))

# Calculate and plot the ALE main effects
par(mfrow = c(2,2))
ALE.1 <- ALEPlot(x2, fit$bart, pred.fun = my_pred_func, J = 1, K = 100)
ALE.2 <- ALEPlot(x2, fit$bart, pred.fun = my_pred_func, J = 2, K = 100)
ALE.3 <- ALEPlot(x2, fit$bart, pred.fun = my_pred_func, J = 3, K = 100)
ALE.4 <- ALEPlot(x2, fit$bart, pred.fun = my_pred_func, J = 4, K = 100)

# Manually plot
par(mfrow = c(2,2))
plot(ALE.1$x.values, ALE.1$f.values, type = "l", xlab = "PM25", ylab = "ALE", xlim = c(10, 20), ylim = c(-0.02, 0.02))
plot(ALE.2$x.values, ALE.2$f.values, type = "l", xlab = "NO2", ylab = "ALE")
plot(ALE.3$x.values, ALE.3$f.values, type = "l", xlab = "O3", ylab = "ALE")
plot(ALE.4$x.values, ALE.4$f.values, type = "l", xlab = "Tmax", ylab = "ALE")


# Calculate and plot the ALE second order effects including temperature
par(mfrow = c(1,3))
ALE.41 <- ALEPlot(x2, fit$bart, pred.fun = my_pred_func, J = c(4,1), K = 100, NA.plot = FALSE)
ALE.42 <- ALEPlot(x2, fit$bart, pred.fun = my_pred_func, J = c(4,2), K = 100, NA.plot = FALSE)
ALE.43 <- ALEPlot(x2, fit$bart, pred.fun = my_pred_func, J = c(4,3), K = 100, NA.plot = FALSE)

# Manually plot
par(mfrow = c(1,3))
fields::image.plot(ALE.41$x.values[[1]], ALE.41$x.values[[2]], ALE.41$f.values, 
                   xlab = "Tmax", ylab = "PM2.5", col = hcl.colors(10, palette = 'Blues'))
contour(ALE.41$x.values[[1]], ALE.41$x.values[[2]], ALE.41$f.values, add = TRUE, drawlabels = TRUE)

fields::image.plot(ALE.42$x.values[[1]], ALE.42$x.values[[2]], ALE.42$f.values, 
                   xlab = "Tmax", ylab = "NO2", col = hcl.colors(10, palette = 'Blues'))
contour(ALE.42$x.values[[1]], ALE.42$x.values[[2]], ALE.42$f.values, add = TRUE, drawlabels = TRUE)

fields::image.plot(ALE.43$x.values[[1]], ALE.43$x.values[[2]], ALE.43$f.values, 
                   xlab = "Tmax", ylab = "O3", col = hcl.colors(10, palette = 'Blues'))
contour(ALE.43$x.values[[1]], ALE.43$x.values[[2]], ALE.43$f.values, add = TRUE, drawlabels = TRUE)

filter(data, Tmax <= 290 & PM25 >= 10) |> nrow()

# Manually plot with main effects included
par(mfrow = c(1,3))
f.41 <- sweep(sweep(ALE.41$f.values, 1, ALE.4$f.values, "+"), 2, ALE.1$f.values, "+")
fields::image.plot(ALE.41$x.values[[1]], ALE.41$x.values[[2]], f.41, 
                   xlab = "Tmax", ylab = "PM2.5", col = hcl.colors(10, palette = 'Blues'))
#contour(ALE.41$x.values[[1]], ALE.41$x.values[[2]], f.41, add = TRUE, drawlabels=TRUE)

f.42 <- sweep(sweep(ALE.42$f.values, 1, ALE.4$f.values, "+"), 2, ALE.2$f.values, "+")
fields::image.plot(ALE.42$x.values[[1]], ALE.42$x.values[[2]], f.42, 
                   xlab = "Tmax", ylab = "NO2", col = hcl.colors(10, palette = 'Blues'))
#contour(ALE.42$x.values[[1]], ALE.42$x.values[[2]], f.42, add = TRUE,drawlabels=TRUE)

f.43 <- sweep(sweep(ALE.43$f.values, 1, ALE.4$f.values, "+"), 2, ALE.3$f.values, "+")
fields::image.plot(ALE.43$x.values[[1]], ALE.43$x.values[[2]], f.43, 
                   xlab = "Tmax", ylab = "O3", col = hcl.colors(10, palette = 'Blues'))
#contour(ALE.43$x.values[[1]], ALE.43$x.values[[2]], f.43, add = TRUE, drawlabels=TRUE)



# Store Output ------------------------------------------------------------
write_rds(fit, paste0(dir, '03. Working Code/02. Analysis/resp_05-07.rds'))


