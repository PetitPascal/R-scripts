#-------------------------------------------------------------------------------
## Reproducible & Generalisable DLNM Script
# - Simulates time-series exposure-response data with lagged effects
# - Fits DLNM
# - Tests assumptions: residual autocorrelation, overdispersion, model fit
# - Visualises 3D and 2D exposure-lag-response surfaces
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup
rm(list = ls())
set.seed(123)

required_pkgs <- c("dlnm", "splines", "MASS", "dplyr", "ggplot2","tidyr", "tseries", "lmtest", "mgcv")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

#--------------------------------------------
## Step 2: Simulating daily time-series data
n_days<-730 # 2 years of daily data

# Time index
day_seq<-seq.Date(as.Date("2020-01-01"), by = "day", length.out = n_days)

# Exposure: daily temperature (Celsius) — seasonal pattern + noise
temp<-15 + 12 * sin(2 * pi * (1:n_days) / 365) + rnorm(n_days, 0, 3)

# Confounders
dow<-as.factor(weekdays(day_seq)) # day of week
trend<-1:n_days / n_days  # long-term trend

# Simulating lagged exposure effect (true: lag 0-7, peak at lag 2)
max_lag<-14
lag_wts<-exp(-0.4 * (0:max_lag)) # decaying lag weights
lag_wts<-lag_wts/sum(lag_wts)

# Building lagged temperature matrix
temp_lagmat<-matrix(NA, nrow = n_days, ncol = max_lag + 1)
for (l in 0:max_lag) {
  temp_lagmat[(l+1):n_days, l+1]<-temp[1:(n_days-l)]
}

# True exposure-response: U-shaped (harm at extremes)
f_true<-function(x) 0.003 * (x - 15)^2

# Weighted lagged effect
h_lag<-rowSums(sweep(temp_lagmat, 2, lag_wts, "*"),na.rm = TRUE)
h_lag_effect<-f_true(h_lag)

# Outcome: daily count (Poisson) — mortality/hospital admissions
log_mu<-3.5 + 0.5 * trend + h_lag_effect

# Day-of-week effect
dow_effect<-c(Mon=0, Tue=0.05, Wed=0.05, Thu=0.03, Fri=0.1, Sat=-0.1, Sun=-0.15)
dow_effect<-c(lun.=0, mar.=0.05, mer.=0.05, jeu.=0.03, ven.=0.1, sam.=-0.1, dim.=-0.15)
log_mu<-log_mu + dow_effect[as.character(format(day_seq, "%a"))]
y_count<-rpois(n_days, exp(log_mu))

ts_data<-data.frame(date=day_seq,
  y=y_count,
  temp=temp,
  trend=trend,
  dow=dow)

# Removing rows with NA due to lag construction
ts_data<-ts_data[(max_lag + 1):n_days, ]
temp_trimmed<-temp[(max_lag + 1):n_days]
cat("Time series length after trimming:", nrow(ts_data), "days\n")

#--------------------------------------------
## Step 3: Assumption / EDA checks

# 3a. Outcome distribution
summary(ts_data$y)
cat("Variance / Mean (expect ~1 for Poisson):",
    round(var(ts_data$y) / mean(ts_data$y), 2), "\n")

# 3b. Autocorrelation in outcome
par(mfrow = c(1, 2))
acf(ts_data$y,  main = "ACF: daily count")
pacf(ts_data$y, main = "PACF: daily count")
par(mfrow = c(1, 1))

# 3c. Exposure distribution and temporal trend
ggplot(ts_data, aes(x = date, y = temp)) +
  geom_line(color = "#2166ac", alpha = 0.7) +
  labs(title = "Daily temperature over time", x = "", y = "Temperature (°C)") +
  theme_bw()

# 3d. Stationarity of exposure (ADF test)
adf_res <- tseries::adf.test(ts_data$temp)
cat("\nADF test (stationarity of exposure):",
    "statistic =", round(adf_res$statistic, 3),
    "| p =", round(adf_res$p.value, 4), "\n")

#--------------------------------------------
## Step 4: Building cross-basis matrix (DLNM)

# Cross-basis combining:
#   - exposure-response function (argvar): natural spline in exposure
#   - lag-response function (arglag): natural spline over lags

cb_temp <- dlnm::crossbasis(ts_data$temp,
  lag=max_lag,
  argvar=list(fun = "ns", df = 4), # 4 df for exposure dimension
  arglag=list(fun = "ns", df = 4)) # 4 df for lag dimension

cat("\nCross-basis dimensions:", dim(cb_temp), "\n")

#--------------------------------------------
## Step 5: Fitting DLNM model

# Poisson GLM with cross-basis + confounders
# Spline for long-term trend (time) to control seasonality

dlnm_fit<-glm(y ~ cb_temp +
    ns(trend, df = 6) +  # flexible time trend
    dow, # day-of-week
  data=ts_data,
  family=quasipoisson()) # quasipoisson to handle overdispersion

summary(dlnm_fit)$coefficients[1:5, ]

#--------------------------------------------
## Step 6: Assumption checks on fitted model

# 6a. Overdispersion (quasipoisson dispersion parameter)
disp_param<-summary(dlnm_fit)$dispersion
cat("\nDispersion parameter (should be ~1 for Poisson):",
    round(disp_param, 3), "\n")
if (disp_param > 1.5)
  cat("Overdispersion present — quasiPoisson family is appropriate.\n")

# 6b. Residual autocorrelation (Durbin-Watson)
dw_test <- lmtest::dwtest(dlnm_fit)
cat("\nDurbin-Watson test: DW =", round(dw_test$statistic, 3),
    "| p =", round(dw_test$p.value, 4),
    "(p < 0.05 = residual autocorrelation present)\n")

# If autocorrelation is severe, add AR(1) term or use ARIMA-type approach
acf(residuals(dlnm_fit, type = "deviance"),
    main = "ACF of deviance residuals (DLNM)")

# 6c. Residuals vs fitted
plot(fitted(dlnm_fit), residuals(dlnm_fit, type = "deviance"),
     xlab = "Fitted", ylab = "Deviance residuals",
     main = "Residuals vs Fitted (DLNM)")
abline(h = 0, lty = 2, col = "red")

#--------------------------------------------
## Step 7: Predicting exposure-lag-response surface

# Reference: median temperature
ref_temp<-median(ts_data$temp)

pred_dlnm <- dlnm::crosspred(basis=cb_temp,
  model=dlnm_fit,
  at=seq(min(ts_data$temp), max(ts_data$temp), length.out = 50),
  cen=ref_temp,
  cumul=TRUE) # also compute cumulative effects

#--------------------------------------------
## Step 8: Visualizing results

# 8a. 3D surface plot
plot(pred_dlnm, xlab = "Temperature (°C)", zlab = "RR",
     main = "DLNM: exposure-lag-response surface",
     theta = 40, phi = 30, col = "lightblue")

# 8b. Overall cumulative exposure-response (summed over all lags)
plot(pred_dlnm, "overall", xlab = "Temperature (°C)", ylab = "Cumulative RR",
     main = "Overall cumulative effect (ref = median temp)",
     col = "#2166ac", lwd = 2, ci.arg = list(col = "#A6DDCE"))

# # 8c. Lag-response at a specific temperature (e.g. 30°C = high)
# plot(pred_dlnm, "slices", var = 30, col = "#d6604d", lwd = 2,
#      main = "Lag-response at 30°C",
#      xlab = "Lag (days)", ylab = "RR")

# 8d. ggplot version of overall cumulative effect
cum_df<-data.frame(temp=pred_dlnm$predvar,
  rr=pred_dlnm$allRRfit,
  lo=pred_dlnm$allRRlow,
  hi=pred_dlnm$allRRhigh)

ggplot(cum_df, aes(x = temp, y = rr, ymin = lo, ymax = hi)) +
  geom_ribbon(fill = "#A6DDCE", alpha = 0.5) +
  geom_line(linewidth = 1.2, color = "#2166ac") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = ref_temp, linetype = "dotted", color = "grey50") +
  labs(title = "DLNM: cumulative exposure-response (ref = median temperature)",
       x = "Temperature (°C)", y = "Cumulative RR (95% CI)") +
  theme_bw(base_size = 14)

#--------------------------------------------
## Step 9: Lag-specific effects at selected temperatures
cold_temp<-quantile(ts_data$temp, 0.05)
hot_temp<-quantile(ts_data$temp, 0.95)

lag_df <- data.frame(lag=0:max_lag,
  rr_cold=pred_dlnm$matRRfit[which.min(abs(pred_dlnm$predvar - cold_temp)), ],
  lo_cold=pred_dlnm$matRRlow[which.min(abs(pred_dlnm$predvar - cold_temp)), ],
  hi_cold=pred_dlnm$matRRhigh[which.min(abs(pred_dlnm$predvar - cold_temp)), ],
  rr_hot=pred_dlnm$matRRfit[which.min(abs(pred_dlnm$predvar - hot_temp)), ],
  lo_hot=pred_dlnm$matRRlow[which.min(abs(pred_dlnm$predvar - hot_temp)), ],
  hi_hot=pred_dlnm$matRRhigh[which.min(abs(pred_dlnm$predvar - hot_temp)), ])

lag_df %>%
  tidyr::pivot_longer(-lag, names_to = c("stat","temp_cat"),
                      names_sep = "_") %>%
  tidyr::pivot_wider(names_from = stat, values_from = value) %>%
  ggplot(aes(x = lag, y = rr, ymin = lo, ymax = hi,
             color = temp_cat, fill = temp_cat)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c(cold = "#4393c3", hot = "#d6604d"),
                     labels = c("Cold (5th pct)", "Hot (95th pct)")) +
  scale_fill_manual(values  = c(cold = "#4393c3", hot = "#d6604d"),
                    labels  = c("Cold (5th pct)", "Hot (95th pct)")) +
  labs(title = "Lag-specific RR at extreme temperatures",
       x = "Lag (days)", y = "RR (95% CI)",
       color = "Temperature", fill = "Temperature") +
  theme_bw(base_size = 14)
