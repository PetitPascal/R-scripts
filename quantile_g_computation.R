#-------------------------------------------------------------------------------
## Reproducible & Generalisable G-Computation Script
# - Simulates data (mixture exposures, covariates, binary/continuous outcome)
# - Implements parametric g-computation manually and via gfoRmula / qgcomp
# - Tests assumptions (positivity, model fit, bootstrap CIs)
# - Includes run_gcomp_pipeline() wrapper
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup
rm(list = ls())
set.seed(123)

required_pkgs<-c("MASS", "dplyr", "tidyr", "ggplot2", "boot","qgcomp", "margins", "broom")

is_installed<-required_pkgs %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed == FALSE)){
  install.packages(required_pkgs[!is_installed],repos = "http://cran.us.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

#--------------------------------------------
## Step 2: Simulating data
n<-600
K<-6
rho<-0.45
Sigma<-rho ^ as.matrix(dist(1:K))
raw<-MASS::mvrnorm(n, mu = rep(0, K), Sigma=Sigma)
expo_df<-as.data.frame(exp(raw / 3))
colnames(expo_df)<-paste0("X", seq_len(K))
expo_names<-colnames(expo_df)

cov_df<-data.frame(age=rnorm(n, 50, 10),sex=rbinom(n, 1, 0.5),bmi=rnorm(n, 26, 4))

true_betas<-c(0.4, 0.3, 0.2, -0.1, 0.05, 0.05)
h_z<-as.numeric(as.matrix(log1p(expo_df)) %*% true_betas)
linpred<--1 + 0.015 * cov_df$age + 0.4 * cov_df$sex + 0.05 * cov_df$bmi + h_z
prob<-plogis(linpred)
y_bin<-rbinom(n, 1, prob)
y_cont<-linpred + rnorm(n, 0, 1.5)

sim_data<-cbind(data.frame(y_bin = y_bin, y_cont = y_cont), cov_df,
                as.data.frame(log1p(expo_df)))  # log-transform exposures

#--------------------------------------------
## Step 3: Assumption checks

# 3a. Positivity: ensuring overlap in exposure distribution across outcome groups

pos_check<-sim_data %>%
  group_by(y_bin) %>%
  summarise(across(all_of(expo_names), mean))
pos_check

# 3b. Correlation among exposures
corr_expo<-cor(sim_data[, expo_names], method = "spearman")
round(max(abs(corr_expo[upper.tri(corr_expo)])), 3)

# 3c. Checking for near-perfect separation (for binary outcome)
tab_check<-table(sim_data$y_bin)
tab_check
if (min(tab_check) < 20) warning("Very few events/non-events — model may not converge reliably.")

#--------------------------------------------
## Step 4: Parametric g-computation (manual)

# 4a. Fitting outcome model (logistic for binary)
form_bin<-as.formula(paste("y_bin ~", paste(c(expo_names, "age", "sex", "bmi"), collapse = " + ")))
fit_bin <-glm(form_bin, data = sim_data, family = binomial())

# Checking model fit
summary(fit_bin)

# Hosmer-Lemeshow goodness-of-fit (for binary)
hl_test<-function(model, g = 10) {
  obs<-model$y
  pred<-fitted(model)
  cuts<-quantile(pred, probs = seq(0, 1, 1/g))
  grp<-cut(pred, breaks = cuts, include.lowest = TRUE)
  tab<-tapply(obs, grp, function(x) c(obs = sum(x), exp = length(x) * mean(x), n = length(x)))
  chi2<-sum(sapply(tab, function(x) (x["obs"] - x["exp"])^2 / (x["exp"] * (1 - x["exp"]/x["n"]))))
  p<-pchisq(chi2, df = g - 2, lower.tail = FALSE)
  cat("\nHosmer-Lemeshow test: chi2 =", round(chi2, 3), ", p =", round(p, 3),
      "(p > 0.05 = adequate fit)\n")
}
hl_test(fit_bin)

# 4b. Counterfactual datasets
q25<-sapply(sim_data[, expo_names], quantile, probs = 0.25)
q75<-sapply(sim_data[, expo_names], quantile, probs = 0.75)

data_low<-sim_data; data_low[, expo_names]<-matrix(q25, nrow = n, ncol = K, byrow = TRUE)
data_high<-sim_data; data_high[, expo_names]<-matrix(q75, nrow = n, ncol = K, byrow = TRUE)

# 4c. Marginal predictions under each scenario
p_low<-mean(predict(fit_bin, newdata = data_low,  type = "response"))
p_high<-mean(predict(fit_bin, newdata = data_high, type = "response"))
RD<-p_high - p_low
RR<-p_high / p_low

cat("\nG-computation results\n")
cat("P(Y=1 | all X at Q25):", round(p_low,  4), "\n")
cat("P(Y=1 | all X at Q75):", round(p_high, 4), "\n")
cat("Risk Difference (RD):", round(RD, 4), "\n")
cat("Risk Ratio (RR):     ", round(RR, 4), "\n")

#--------------------------------------------
## Step 5: Bootstrap confidence intervals
gcomp_boot<-function(data, indices, expo_names, outcome_col,covariate_cols, K) {
  d<-data[indices, ]
  form<-as.formula(paste(outcome_col, "~",paste(c(expo_names, covariate_cols), collapse = " + ")))
  fit<-glm(form, data = d, family = binomial())
  q25<-sapply(d[, expo_names], quantile, probs = 0.25)
  q75<-sapply(d[, expo_names], quantile, probs = 0.75)
  d_lo<-d; d_lo[, expo_names]<-matrix(q25, nrow = nrow(d), ncol = K, byrow = TRUE)
  d_hi<-d; d_hi[, expo_names]<-matrix(q75, nrow = nrow(d), ncol = K, byrow = TRUE)
  c(RD = mean(predict(fit, d_hi, type = "response")) - mean(predict(fit, d_lo, type = "response")),
    RR = mean(predict(fit, d_hi, type = "response")) / mean(predict(fit, d_lo, type = "response")))
}

set.seed(42)
boot_res<-boot::boot(data = sim_data, statistic = gcomp_boot,
                       R = 500,   # increase to 1000-2000 for publication
                       expo_names = expo_names, outcome_col = "y_bin",
                       covariate_cols = c("age", "sex", "bmi"), K = K)

ci_RD<-boot::boot.ci(boot_res, index = 1, type = "perc")
ci_RR<-boot::boot.ci(boot_res, index = 2, type = "perc")

cat("\nBootstrap 95% CIs\n")
cat("RD:", round(RD, 4), "  95% CI:", round(ci_RD$percent[4], 4), "to", round(ci_RD$percent[5], 4), "\n")
cat("RR:", round(RR, 4), "  95% CI:", round(ci_RR$percent[4], 4), "to", round(ci_RR$percent[5], 4), "\n")

#--------------------------------------------
## Step 6: Exposure-specific effects via qgcomp
form_qgc<-as.formula(paste("y_bin ~", paste(c(expo_names, "age", "sex", "bmi"), collapse = " + ")))
qgc_fit<-qgcomp::qgcomp.noboot(f=form_qgc,
  expnms=expo_names,
  data=sim_data,
  family=binomial(),
  q=4)
qgc_fit
plot(qgc_fit)

# With bootstrap (for valid CIs)
qgc_boot<-qgcomp::qgcomp.boot(f=form_qgc,
  expnms=expo_names,
  data=sim_data,
  family=binomial(),
  q=4,
  B=500,
  seed=123)
qgc_boot
plot(qgc_boot)

#--------------------------------------------
## Step 7: Continuous outcome g-computation
form_cont<-as.formula(paste("y_cont ~", paste(c(expo_names, "age", "sex", "bmi"), collapse = " + ")))
fit_cont <-lm(form_cont, data = sim_data)

# Residual checks
par(mfrow = c(2, 2)); plot(fit_cont); par(mfrow = c(1, 1))

# Checking normality of residuals
sw_test<-shapiro.test(residuals(fit_cont))
sw_test

# Homoscedasticity (Breusch-Pagan)
if (requireNamespace("lmtest", quietly = TRUE) && requireNamespace("car", quietly = TRUE)) {
  bp_test<-lmtest::bptest(fit_cont)
  bp_test
}

mu_low <-mean(predict(fit_cont, newdata = data_low))
mu_high<-mean(predict(fit_cont, newdata = data_high))
cat("\nG-computation (continuous)\n")
cat("E[Y | X at Q25]:", round(mu_low, 4), "\n")
cat("E[Y | X at Q75]:", round(mu_high, 4), "\n")
cat("Mean difference:", round(mu_high - mu_low, 4), "\n")

#--------------------------------------------
## Step 8: Visualisation

# Bootstrap distribution of RD
boot_rd_df<-data.frame(RD = boot_res$t[, 1])
ggplot(boot_rd_df, aes(x = RD)) +
  geom_histogram(fill = "#A6DDCE", color = "black", bins = 40) +
  geom_vline(xintercept = RD, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = c(ci_RD$percent[4], ci_RD$percent[5]),
             color = "navy", linetype = "dotted") +
  labs(title = "Bootstrap distribution: Risk Difference (g-computation)",
       x = "Risk Difference", y = "Count") +
  theme_bw(base_size = 14)

# Marginal predictions plot across quantile scenarios
q_seq<-seq(0.10, 0.90, by = 0.10)
pred_seq<-sapply(q_seq, function(q) {
  q_val<-sapply(sim_data[, expo_names], quantile, probs = q)
  d_q  <-sim_data
  d_q[, expo_names]<-matrix(q_val, nrow = n, ncol = K, byrow = TRUE)
  mean(predict(fit_bin, newdata = d_q, type = "response"))
})

pred_df<-data.frame(quantile = q_seq, predicted_prob = pred_seq)
ggplot(pred_df, aes(x = quantile, y = predicted_prob)) +
  geom_line(linewidth = 1.2, color = "#2166ac") +
  geom_point(size = 3, color = "#2166ac") +
  scale_x_continuous(labels = scales::percent) +
  labs(title = "G-computation: predicted P(Y=1) across joint exposure quantiles",
       x = "Quantile of all exposures (joint shift)", y = "Predicted probability") +
  theme_bw(base_size = 14)

#--------------------------------------------
## Step 9: Reusable pipeline
run_gcomp_pipeline<-function(data, outcome_col, expo_cols, covariate_cols = NULL,
                               family = c("binomial", "gaussian"),
                               q_low = 0.25, q_high = 0.75,
                               n_boot = 500, seed = 2025,
                               log_transform = TRUE) {
  family<-match.arg(family)
  set.seed(seed)
  df<-data
  
  if (log_transform) df[, expo_cols]<-log1p(df[, expo_cols])
  
  all_preds<-c(expo_cols, covariate_cols)
  form<-as.formula(paste(outcome_col, "~", paste(all_preds, collapse = " + ")))
  fam <-if (family == "binomial") binomial() else gaussian()
  fit <-glm(form, data = df, family = fam)
  
  K <-length(expo_cols)
  nn<-nrow(df)
  q25<-sapply(df[, expo_cols], quantile, probs = q_low)
  q75<-sapply(df[, expo_cols], quantile, probs = q_high)
  
  d_lo<-df; d_lo[, expo_cols]<-matrix(q25, nrow = nn, ncol = K, byrow = TRUE)
  d_hi<-df; d_hi[, expo_cols]<-matrix(q75, nrow = nn, ncol = K, byrow = TRUE)
  
  if (family == "binomial") {
    p_lo<-mean(predict(fit, d_lo, type = "response"))
    p_hi<-mean(predict(fit, d_hi, type = "response"))
    est <-c(p_low = p_lo, p_high = p_hi, RD = p_hi - p_lo, RR = p_hi / p_lo)
  } else {
    m_lo<-mean(predict(fit, d_lo))
    m_hi<-mean(predict(fit, d_hi))
    est <-c(mean_low = m_lo, mean_high = m_hi, diff = m_hi - m_lo)
  }
  
  # Bootstrap
  boot_fn<-function(d, idx) {
    db <-d[idx, ]
    fb <-tryCatch(glm(form, data = db, family = fam), error = function(e) NULL)
    if (is.null(fb)) return(rep(NA, length(est)))
    q25b<-sapply(db[, expo_cols], quantile, probs = q_low)
    q75b<-sapply(db[, expo_cols], quantile, probs = q_high)
    dlob<-db; dlob[, expo_cols]<-matrix(q25b, nrow = nrow(db), ncol = K, byrow = TRUE)
    dhib<-db; dhib[, expo_cols]<-matrix(q75b, nrow = nrow(db), ncol = K, byrow = TRUE)
    if (family == "binomial") {
      p_lo_b<-mean(predict(fb, dlob, type = "response"))
      p_hi_b<-mean(predict(fb, dhib, type = "response"))
      return(c(p_lo_b, p_hi_b, p_hi_b - p_lo_b, p_hi_b / p_lo_b))
    } else {
      ml<-mean(predict(fb, dlob)); mh<-mean(predict(fb, dhib))
      return(c(ml, mh, mh - ml))
    }
  }
  
  br  <-boot::boot(df, boot_fn, R = n_boot)
  ci_l<-apply(br$t, 2, quantile, 0.025, na.rm = TRUE)
  ci_u<-apply(br$t, 2, quantile, 0.975, na.rm = TRUE)
  
  result_df<-data.frame(estimate = est, ci_low = ci_l, ci_high = ci_u)
  cat("\nG-computation results\n")
  print(round(result_df, 4))
  return(list(fit = fit, estimates = result_df, boot = br))
}

res<-run_gcomp_pipeline(
  data          = sim_data,
  outcome_col   = "y_bin",
  expo_cols     = expo_names,
  covariate_cols = c("age", "sex", "bmi"),
  family        = "binomial",
  n_boot        = 500,
  log_transform = FALSE) # already log-transformed above

res
