#-------------------------------------------------------------------------------
## Reproducible & Generalisable GAMM Script
# - Simulates longitudinal data with non-linear exposure effects
# - Fits GAMMs with random intercepts
# - Tests assumptions: residual normality, homoscedasticity, concurvity
# - Produces partial effect plots, saves results
# - Includes run_gamm_pipeline() wrapper
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup
rm(list = ls())
set.seed(123)

required_pkgs<-c("mgcv", "MASS", "dplyr", "ggplot2", "tidyr","gratia", "lmtest", "performance")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

#--------------------------------------------
## Step 2: Simulating longitudinal data
n_id<-120
n_time<-4
ID<-rep(seq_len(n_id), each = n_time)
time<-rep(0:(n_time-1), times = n_id)

# Random subject intercepts
u_i<-rnorm(n_id, 0, 1.5)[ID]

# Non-linear exposure effect
exposure<-rnorm(n_id * n_time, 5, 2)
h_expo<-2 * sin(exposure / 3)   # true non-linear function

# Covariate
age<-rnorm(n_id, 50, 10)[ID]
sex<-rbinom(n_id, 1, 0.5)[ID]

# Continuous outcome
y<-5 + h_expo + 0.05 * age + 0.5 * sex + u_i + rnorm(n_id * n_time, 0, 1)

sim_long<-data.frame(id=factor(ID),
  time=time,
  exposure=exposure,
  age=age,
  sex=sex,
  y=y)

#--------------------------------------------
## Step 3: Assumption checks

# 3a. Outcome distribution
hist(sim_long$y, main = "Outcome distribution", xlab = "y", col = "#A6DDCE")

# 3b. Visualizing exposure vs. outcome (hints at non-linearity)
ggplot(sim_long, aes(x = exposure, y = y)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  geom_smooth(method = "lm",    se = TRUE, color = "blue", linetype = "dashed") +
  labs(title = "Exposure vs outcome: LOESS vs linear", x = "Exposure", y = "y") +
  theme_bw()

# 3c. Checking for sufficient variation in exposure
cat("Exposure summary:\n"); print(summary(sim_long$exposure))
cat("SD of exposure:", round(sd(sim_long$exposure), 3), "\n")

# 3d. Number of observations per subject
obs_per_id<-table(sim_long$id)
cat("Observations per subject: min =", min(obs_per_id),
    ", max =", max(obs_per_id), "\n")

#--------------------------------------------
## Step 4: Fitting GAMM
# Using mgcv::bam for large data or gamm4::gamm4 for lme4-based random effects
# Here: mgcv::gam with bs="re" for random intercept

gamm_fit<-mgcv::gam(
  y ~ s(exposure, bs = "cr", k = 8) +   # smooth for exposure (cubic regression spline)
    s(time,     bs = "cr", k = 4) +   # smooth for time
    age + sex +
    s(id, bs = "re"),                  # random intercept for subject
  data   = sim_long,
  method = "REML"   # REML for smoothing parameter selection
)

summary(gamm_fit)
broom::tidy(gamm_fit)

#--------------------------------------------
## Step 5: Assumption testing

# 5a. Residual diagnostics
par(mfrow = c(2, 2))
mgcv::gam.check(gamm_fit)
par(mfrow = c(1, 1))

# 5b. Normality of residuals
shapiro.test(residuals(gamm_fit, type = "deviance"))

# 5c. Concurvity (analogous to multicollinearity for GAMs)
conc<-mgcv::concurvity(gamm_fit, full = FALSE)
round(conc$worst, 3)

# 5d. Basis dimension adequacy (k-index; should be > 1)
mgcv::k.check(gamm_fit)

# 5e. Testing non-linearity: comparing with linear model
lm_fit<-lm(y ~ exposure + time + age + sex, data = sim_long)
anova(lm_fit, gamm_fit, test = "F")

#--------------------------------------------
## Step 6: Visualizing smooth effects

# Using gratia
gratia::draw(gamm_fit, select = 1:2)   # exposure and time smooths

# Manual ggplot
smooth_df<-gratia::smooth_estimates(gamm_fit, smooth = "s(exposure)")
ggplot(smooth_df, aes(x = exposure, y = est)) +
  geom_ribbon(aes(ymin = est - 2*se, ymax = est + 2*se), alpha = 0.25, fill = "#A6DDCE") +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "GAMM: Partial effect of exposure",
       x = "Exposure", y = "s(exposure)") +
  theme_bw(base_size = 14)

#--------------------------------------------
## Step 7: Predict and plot fitted values-
sim_long$fitted_gamm<-fitted(gamm_fit)
ggplot(sim_long, aes(x = fitted_gamm, y = y)) +
  geom_point(alpha = 0.3, color = "#2166ac") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Observed vs fitted values (GAMM)",
       x = "Fitted", y = "Observed") +
  theme_bw(base_size = 14)

#--------------------------------------------
## Step 8: Binary outcome GAMM (logistic)
y_bin<-rbinom(nrow(sim_long), 1, plogis(sim_long$y / 5))
sim_long$y_bin<-y_bin

gamm_bin<-mgcv::gam(
  y_bin ~ s(exposure, bs = "cr", k = 8) +
    s(time, bs = "cr", k = 4) +
    age + sex +
    s(id, bs = "re"),
  data=sim_long,
  family=binomial(),
  method="REML")

summary(gamm_bin)
broom::tidy(gamm_bin)
mgcv::gam.check(gamm_bin)

#--------------------------------------------
## Step 9: Reusable pipeline
run_gamm_pipeline<-function(data, outcome_col, smooth_vars, linear_vars = NULL,
                              id_col = NULL, k =8, family = gaussian(),
                              method = "REML") {
  # Building smooth terms
  smooth_terms<-paste0("s(", smooth_vars, ", bs='cr', k=", k, ")", collapse = " + ")
  
  # Building linear terms
  linear_terms<-if (!is.null(linear_vars)) paste(linear_vars, collapse = " + ") else ""
  
  # Random intercept
  re_term<-if (!is.null(id_col)) paste0("s(", id_col, ", bs='re')") else ""
  
  all_terms<-paste(c(smooth_terms, linear_terms, re_term)[nchar(c(smooth_terms, linear_terms, re_term)) > 0],
                     collapse = " + ")
  form<-as.formula(paste(outcome_col, "~", all_terms))
  
  fit<-mgcv::gam(form, data = data, family = family, method = method)
  
  # Diagnostics
  cat("\n--- GAMM Summary ---\n")
  print(summary(fit))
  cat("\n--- Concurvity ---\n")
  print(round(mgcv::concurvity(fit, full = FALSE)$worst, 3))
  cat("\n--- k-index ---\n")
  mgcv::k.check(fit)
  
  return(fit)
}

fit_pipe<-run_gamm_pipeline(data=sim_long,
  outcome_col="y",
  smooth_vars=c("exposure", "time"),
  linear_vars=c("age", "sex"),
  id_col="id",
  k=4,
  family=gaussian())
