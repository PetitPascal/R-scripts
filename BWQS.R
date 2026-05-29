#-------------------------------------------------------------------------------
## Reproducible & Generalisable Bayesian Weighted Quantile Sum (BWQS) Script
# - Simulates correlated mixture exposures, covariates, and a binary outcome
# - Fits BWQS model using the bwqs package
# - Tests assumptions, extracts weights, plots results
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup
rm(list = ls())
set.seed(123)

required_pkgs<-c("BWQS", "MASS", "dplyr", "ggplot2", "corrplot","tidyr", "forcats", "rstan")

if(!("BWQS"%in%.packages(all.available=TRUE))){
  devtools::install_github("ElenaColicino/bwqs", build_vignettes = TRUE)
}

is_installed<-required_pkgs %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed == FALSE)){
  install.packages(required_pkgs[!is_installed])
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

#--------------------------------------------
## Step 2: Simulating data
n<-400
K<-8   # number of mixture components

rho<-0.5
Sigma<-rho ^ as.matrix(dist(1:K))
raw<-MASS::mvrnorm(n, mu = rep(0, K), Sigma = Sigma)

expo_df<-as.data.frame(exp(raw / 3))
colnames(expo_df)<-paste0("chem_", seq_len(K))

cov_df<-data.frame(age=rnorm(n, 50, 12),sex=rbinom(n, 1, 0.5))

true_w<-c(0.35, 0.25, 0.15, 0.10, rep(0.0375, 4))  # sums to 1
h_z<-as.numeric(as.matrix(log1p(expo_df)) %*% true_w)
linpred<- -1.5 + 0.02 * cov_df$age + 0.3 * cov_df$sex + 2 * h_z
prob<-plogis(linpred)
y<-rbinom(n, 1, prob)

sim_data<-cbind(data.frame(y = y), cov_df, expo_df)
expo_names<-paste0("chem_", seq_len(K))

#--------------------------------------------
## Step 3: EDA — correlation among exposures
corr_mat<-cor(log1p(sim_data[, expo_names]), method = "spearman")
corrplot::corrplot.mixed(corr_mat, tl.cex = 0.8, number.cex = 0.7,
                         main = "Spearman correlation: mixture components")

# Checking exposure distributions
par(mfrow = c(2, 4))
for (v in expo_names) hist(log1p(sim_data[[v]]), main = v, xlab = "log1p")
par(mfrow = c(1, 1))

#--------------------------------------------
## Step 4: Quantizing exposures (required by BWQS)
# bwqs expects integer quantile-ranked exposures (1 to q)
quantize<-function(x, q = 4) as.integer(cut(x, breaks = quantile(x, probs = seq(0, 1, 1/q),
                                                                 na.rm = TRUE),
                                            include.lowest = TRUE))
q_val<-4
expo_q<-as.data.frame(lapply(sim_data[, expo_names], quantize, q = q_val))
colnames(expo_q)<-expo_names

#--------------------------------------------
## Step 5: Fitting BWQS model
# bwqs() arguments:
#   formula : outcome ~ covariates (mixture handled separately via mix_name)
#   mix_name: character vector of mixture variable names
#   data: data containing quantized exposures + covariates + outcome
#   q: number of quantiles already applied
#   chains/iter: Stan MCMC settings — increase for final analysis
#   family: "binomial" or "gaussian"

bwqs_data<-cbind(data.frame(y = y), cov_df, expo_q)

bwqs_fit<-BWQS::bwqs(formula= y ~ age + sex,
                     mix_name=expo_names,
                     data=bwqs_data,
                     q=NULL,        # already quantized above
                     family="binomial",
                     chains=2,
                     c_int = c(0.025, 0.975),
                     prior = "None",
                     thin = 3,
                     Dalp = NULL,
                     start_value = NULL,
                     iter=2000,        # increase to 4000+ for publication
                     seed=123)

#--------------------------------------------
## Step 6: Summarizing posterior
bwqs_fit

# Extracting the weight
Weight<-bwqs_fit$summary_fit[grep("^W_", rownames(bwqs_fit$summary_fit)), , drop = FALSE]
Weight<- as_tibble(data.frame(exposure=rownames(Weight),Weight))
Weight

# Extracting mixture effect
beta_summary<-bwqs_fit$summary_fit[grep("beta_wqs|b_wqs", rownames(bwqs_fit$summary_fit)), ]
beta_summary

#--------------------------------------------
## Step 7: Assumption / diagnostics checks ----

# 7a. MCMC convergence: R-hat < 1.01 for all parameters
rhat_vals<-bwqs_fit$summary_fit[, "Rhat"]
cat("\n--- R-hat diagnostics (should all be < 1.01) ---\n")
print(summary(rhat_vals))
if (any(rhat_vals > 1.01, na.rm = TRUE)) {
  warning("Some R-hat > 1.01: consider more iterations or checking priors.")
}

# 7b. Effective sample size (n_eff)
neff_vals<-bwqs_fit$summary_fit[, "n_eff"]
cat("\n--- Effective sample sizes (n_eff) ---\n")
print(summary(neff_vals))
if (any(neff_vals < 100, na.rm = TRUE)) {
  warning("Some n_eff < 100: poor mixing. Increase iterations.")
}

#--------------------------------------------
## Step 8: Plotting weight
ggplot(Weight, aes(x = reorder(exposure, mean), y = mean)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = NULL, y = "Posterior mean weight") +
  theme_minimal()

ggplot(Weight, aes(x = reorder(exposure, mean), y = mean)) +
  geom_point(col = "steelblue") +
  coord_flip() +
  geom_errorbar(aes(ymin =X2.5., ymax =X97.5.), width = 0.2) +
  coord_flip() +
  labs(x = NULL, y = "Posterior mean weight") +
  theme_minimal()

#--------------------------------------------
## Step 9: Plotting mixture effect
library(bayesplot)

beta1_draws <- as.data.frame(rstan::extract(bwqs_fit$fit)$beta1)
colnames(beta1_draws)<-"beta1"

ggplot(beta1_draws, aes(x = beta1)) +
  geom_density(fill = "tomato", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_minimal() +
  labs(x = expression(beta[1]), y = "Posterior density")
