#-------------------------------------------------------------------------------
## Reproducible & generalizable Bayesian kernel machine regression (BKMR) analysis script
# - simulates a correlated exposure mixture, covariates, and a binary outcome
# - runs BKMR (kmbayes) with variable selection, knots, and posterior summaries
# - produces diagnostics and plots
# - use it as a template for your own dataset by replacing the simulated data
#-------------------------------------------------------------------------------

## -------------------------
## Step 1 - Setup: packages + helpers
## -------------------------

required_pkgs <- c("bkmr", "fields", "openxlsx", "ggplot2", "corrplot", "MASS",
                   "dplyr", "tidyr", "scales", "gridExtra")
for(p in required_pkgs){
  if(!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

set.seed(2025)  # reproducibility seed (change as needed)

## Small helper to standardize a matrix (center & scale)
standardize_mat <- function(m){
  m <- as.matrix(m)
  apply(m, 2, function(x) if(sd(x, na.rm = TRUE) > 0) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) else x)
}

## -------------------------
## Step 2 - Simulating a dataset (replace this block with code to load your own data)
## -------------------------
n <- 500
p_exposures <- 10

# Creating correlated exposures using multivariate normal, then exponentiating to resemble concentrations
Sigma <- 0.6 ^ as.matrix(dist(1:p_exposures))   # exponential decay correlation
raw_expos <- MASS::mvrnorm(n = n, mu = rep(0, p_exposures), Sigma = Sigma)

# Making exposures positive and on a scale like concentrations
expo_df <- as.data.frame(exp(raw_expos / 3))    # scale down variation

# Naming exposures with a prefix 'X_'
colnames(expo_df) <- paste0("X_", seq_len(p_exposures))

# Covariates: continuous age, BMI, and a binary sex variable
cov_df <- data.frame(
  age_years = round(rnorm(n, mean = 50, sd = 12)),
  bmi = rnorm(n, mean = 27, sd = 5),
  sex_male = rbinom(n, 1, p = 0.48)
)

# Creating an outcome (binary) where a non-linear combination of exposures influences probability:
# true underlying function for simulation: h(Z) = sum(w_j * log(X_j+1)) + interaction between X1 and X2
weights <- c(0.4, 0.25, rep(0.05, p_exposures-2))
hZ <- as.numeric(as.matrix(log(expo_df + 1)) %*% weights)
hZ <- hZ + 0.6 * (log(expo_df$X_1 + 1) * log(expo_df$X_2 + 1))  # some interaction

# Linear predictor includes covariate effects
linpred <- -2 + 0.02 * cov_df$age_years + 0.06 * cov_df$bmi + 0.4 * cov_df$sex_male + hZ
prob <- 1 / (1 + exp(-linpred))
outcome_y <- rbinom(n, 1, prob)

# Combining into one data.frame
sim_data <- bind_cols(data.frame(id = seq_len(n), outcome = outcome_y), cov_df, expo_df)

## -------------------------
## Step 3 - Quick EDA: correlation matrix of exposures + pairs plot (small)
## -------------------------

expo_mat <- as.matrix(sim_data %>% select(starts_with("X_")))
corr_mat <- cor(expo_mat, method = "spearman")
corrplot::corrplot.mixed(corr_mat, number.cex = 0.7, tl.cex = 0.8, main = "Spearman correlations: exposures")

## optional pairs: (comment out if too large)
# pairs(log(expo_mat + 1), main = "Pairwise log-exposure relationships")

## -------------------------
## Step 4 - Preparing matrices for BKMR (standardized exposures)
## -------------------------
expo_std <- standardize_mat(log(expo_mat + 1))  # log + standardization
covs_matrix <- as.matrix(sim_data %>% select(age_years, bmi, sex_male))
outcome_vector <- sim_data$outcome

# Names for exposures for plotting labels
expo_names <- colnames(expo_std)

## -------------------------
## Step 5 - Determining knots (cover.design) and groups (optional)
#  
#  knots: choose number based on n and complexity (rule of thumb: 25-100)
## -------------------------

# Choosing 50 knots by default (change as needed)
num_knots <- 50
knots_coords <- fields::cover.design(expo_std, nd = num_knots)$design

## Demonstrating a possible group structure (example: group some exposures)
# you can create a vector length p_exposures specifying group ids or leave NULL
group_vector <- c(1,1,2,3,3,4,4,5,6,7)  # arbitrary example; change if you have domain knowledge

## -------------------------
## Step 6 - Fitting BKMR model
## -------------------------
## NOTES:
# - family = "binomial" for binary outcome
# - iter: set to 5000 for demonstration; consider 10000-30000 for final inference
# - varsel = TRUE enables variable selection (PIPs)
# - est.h = TRUE will estimate h(z) for predictor-response plots
# - control.params can be tuned; r.jump2 controls MCMC mixing for h

set.seed(1234)
bkmr_fit <- kmbayes(
  y = outcome_vector,
  Z = expo_std,
  X = covs_matrix,
  iter = 5000,
  family = "binomial",
  varsel = TRUE,
  knots = knots_coords,
  groups = group_vector,
  est.h = TRUE,
  verbose = TRUE,
  control.params = list(r.jump2 = 0.5)
)

## -------------------------
## Step 7 - Posterior inclusion probabilities (PIPs) and variable selection
## -------------------------
pip_df <- ExtractPIPs(bkmr_fit)
print("Posterior inclusion probabilities (PIPs):")
print(pip_df)

## Visualizing PIPs as a bar plot

# Computing overall variable-level PIP
pip_df$PIP <- pip_df$groupPIP * pip_df$condPIP

# Preparing data for plotting
pip_plot_df <- pip_df[, c("variable", "PIP")]

ggplot(pip_plot_df, aes(x = reorder(variable, PIP), y = PIP)) +
  geom_col(fill="#A6DDCE",col="black") + 
  geom_text(aes(x = reorder(variable, PIP), y = PIP+0.05,label=as.character(signif(PIP,3))),size=5)+
  coord_flip() + 
  ylab("Posterior Inclusion Probability") + xlab("")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0),
        strip.background = element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y = element_text(size=16,color="black"),
        axis.text.x = element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.position="bottom",
        axis.line = element_line(color = "black",size = 0.1, linetype = "solid"))

## -------------------------
## Step 8 - Univariate predictor-response functions: effect of each exposure holding others at median
## -------------------------
pred_univar <- PredictorResponseUnivar(fit = bkmr_fit)

# pred_univar has columns: variable, z, est, se
ggplot(pred_univar, aes(z, est, ymin = est - 1.96 * se, ymax = est + 1.96 * se)) +
  geom_ribbon(alpha = 0.18) +
  geom_line() +
  facet_wrap(~ variable, scales = "free_x") +
  labs(x = "Exposure (standardized)", y = "Estimated h(z)", title = "Univariate predictor-response estimates") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0),
        strip.background = element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y = element_text(size=16,color="black"),
        axis.text.x = element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.position="bottom",
        axis.line = element_line(color = "black",size = 0.1, linetype = "solid"))

## -------------------------
## Step 9 - Overall mixture risk summary
## -------------------------
qs_seq <- seq(0.25, 0.75, by = 0.05)
overall_risks <- OverallRiskSummaries(fit = bkmr_fit, qs = qs_seq, q.fixed = 0.5)

ggplot(overall_risks, aes(quantile, est, ymin = est - 1.96 * sd, ymax = est + 1.96 * sd)) +
  geom_pointrange() + geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Overall mixture risk: change from median (q.fixed=0.5)", x = "Quantile of joint change", y = "Estimated change in h(z)") +
  theme_minimal()+
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0),
        strip.background = element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y = element_text(size=16,color="black"),
        axis.text.x = element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.position="bottom",
        axis.line = element_line(color = "black",size = 0.1, linetype = "solid"))

## -------------------------
## Step 10 - Single-exposure risk summaries (pairwise comparisons)
## -------------------------
single_risks <- SingVarRiskSummaries(
  fit = bkmr_fit,
  y = outcome_vector,
  Z = expo_std,
  X = covs_matrix,
  qs.diff = c(0.25, 0.75),        # compare these quantiles
  q.fixed = c(0.25, 0.50, 0.75)   # settings to fix other exposures
)

# Ploting single-exposure risks
ggplot(single_risks, aes(variable, est, ymin = est - 1.96 * sd, ymax = est + 1.96 * sd, color = factor(q.fixed))) +
  geom_pointrange(position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Single-exposure risk summaries", color = "q.fixed")

## -------------------------
## Step 11 - Diagnostic checks & posterior summaries
## -------------------------
# Basic trace plot for kernel parameters if available
if(!is.null(bkmr_fit$chain$h)){
  matplot(bkmr_fit$chain$h, type = "l", lty = 1, main = "Trace plot: h chain (subset)"); legend("topright", legend = "h", bty = "n")
}

# Saving the fitted object for later use
saveRDS(bkmr_fit, file = "bkmr_fit_example.rds")
message("Saved fitted BKMR object to: bkmr_fit_example.rds")

## -------------------------
## Step 12 - Utility: function to run full analysis on user data
## -------------------------
run_bkmr_pipeline <- function(df, outcome_col, exposure_cols, covariate_cols = NULL,
                              iter = 5000, knots_nd = 50, groups = NULL, seed = 1234){
  set.seed(seed)
  Z <- as.matrix(df[, exposure_cols])
  # applying log+standardization if exposures are positive; if not, just standardize
  Z_trans <- standardize_mat(log(Z + 1))
  X_mat <- if(!is.null(covariate_cols)) as.matrix(df[, covariate_cols]) else NULL
  y_vec <- df[[outcome_col]]
  knots_mat <- fields::cover.design(Z_trans, nd = knots_nd)$design
  fit <- kmbayes(y = y_vec, Z = Z_trans, X = X_mat, iter = iter, family = "binomial",
                 varsel = TRUE, knots = knots_mat, groups = groups, est.h = TRUE, verbose = FALSE)
  return(list(fit = fit, Z_trans = Z_trans))
}

# Example
pipeline_res <- run_bkmr_pipeline(sim_data, outcome_col = "outcome",
                                  exposure_cols = paste0("X_", 1:p_exposures),
                                  covariate_cols = c("age_years", "bmi", "sex_male"),
                                  iter = 5000, knots_nd = 50, groups = group_vector)

## -------------------------
## Notes & recommendations
## -------------------------
cat("Notes:
- This script simulates data for demonstration. To use your own data, replace the simulation block
  and set 'expo_std', 'covs_matrix', and 'outcome_vector' accordingly.
- For final inference increase 'iter' (e.g., 10000-30000) and check mixing (trace plots) and effective
  sample sizes. MCMC heavy computations can take many minutes to hours depending on data and iter.
- Consider sensitivity analyses: different knot numbers, different priors, or excluding influential observations.
- If you have domain knowledge for grouping exposures (e.g., chemical classes), use the 'groups' argument.")
