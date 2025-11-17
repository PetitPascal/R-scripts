#-------------------------------------------------------------------------------
## Restricted Cubic Spline (RCS) pipeline — reproducible & generalisable
# - Simulates data (correlated exposures, covariates, outcome)
# - Fits RCS using rms::lrm (binary) or rms::ols (continuous)
# - Produces spline plot with CIs, tests for nonlinearity, saves outputs
# - Includes run_rcs_pipeline() to apply to your own dataset
#-------------------------------------------------------------------------------

## -------------------------
## Step 1 - Setup
## -------------------------
rm(list = ls())
set.seed(2025)

required_pkgs <- c("MASS", "rms", "ggplot2", "dplyr", "tidyr", "cowplot")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# Creating a helper for safe file paths
safe_filepath <- function(fname) file.path(getwd(), fname)

## -------------------------
## Step 2 - Simulating example data
## -------------------------
n <- 800
p_expo <- 6

# correlated exposures via MVN
rho <- 0.5
Sigma <- rho ^ as.matrix(dist(1:p_expo))
raw <- MASS::mvrnorm(n, mu = rep(0, p_expo), Sigma = Sigma)

exp_df <- as.data.frame(exp(raw / 3))      # making them positive-ish like concentrations
expo_names <- paste0("env_", LETTERS[1:p_expo])
colnames(exp_df) <- expo_names

# covariates
cov_df <- data.frame(
  age_yrs = round(rnorm(n, 52, 12)),
  bmi_kgm2 = rnorm(n, 27, 4),
  male_flag = rbinom(n, 1, 0.48)
)

# creating non-linear true effect for first exposure
h_true <- 0.6 * (log1p(exp_df[[1]]) - 0.5)^2 + 0.35 * log1p(exp_df[[2]])
linpred <- -2 + 0.02 * cov_df$age_yrs + 0.05 * cov_df$bmi_kgm2 + 0.5 * cov_df$male_flag + h_true
prob <- 1 / (1 + exp(-linpred))
y_bin <- rbinom(n, 1, prob)

# continuous outcome example
y_cont <- 10 + 0.3 * cov_df$age_yrs - 0.2 * cov_df$bmi_kgm2 + 2 * (log1p(exp_df[[1]]) ^ 1.5) + rnorm(n, 0, 3)

example_df <- bind_cols(
  tibble::tibble(id = seq_len(n)),
  tibble::tibble(binary_outcome = y_bin, continuous_outcome = y_cont),
  cov_df,
  exp_df
)

head(example_df)

## -------------------------
## Step 3 - Utility: pretty plot theme
## -------------------------
theme_rcs <- function() theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

## -------------------------
## Step 4 - Core function: running one RCS model and plot
## -------------------------
# Parameters:
# - data: data.frame
# - outcome: name of outcome column (string)
# - exposure: name of exposure column to model with rcs (string)
# - covariates: character vector of covariate names
# - family: "binomial" or "gaussian"
# - nk: number of knots for rcs (integer >= 3). If numeric vector provided, used as knots percentiles.
# - ref_pct: percentile (0-1) to pick reference point (OR = 1); if NULL uses mean
# - nsim: grid points to draw curve
# - out_prefix: prefix for saved files

run_rcs_one <- function(data, outcome, exposure, covariates = NULL,
                        family = c("binomial", "gaussian"),
                        nk = 4, ref_pct = 0.25, nsim = 200) {

  family <- match.arg(family)
  
  df <- data %>% select(all_of(c(outcome, exposure, covariates))) %>% na.omit()
  
  dd <- datadist(df)
  options(datadist="dd")
  
  covpart <- if (!is.null(covariates) && length(covariates)>0) paste(covariates, collapse=" + ") else ""
  form_str <- if (covpart=="") {
    paste0(outcome, " ~ rcs(", exposure, ", ", nk, ")")
  } else {
    paste0(outcome, " ~ rcs(", exposure, ", ", nk, ") + ", covpart)
  }
  formula_rcs <- as.formula(form_str)
  
  fit <- if(family=="binomial") lrm(formula_rcs, data=df, x=TRUE, y=TRUE) else ols(formula_rcs, data=df, x=TRUE, y=TRUE)
  
  pred_obj <- if(family=="binomial") {
    do.call(Predict, list(fit, as.name(exposure), fun=exp, type="predictions", np=nsim))
  } else {
    do.call(Predict, list(fit, as.name(exposure), type="predictions", np=nsim))
  }
  
  pred_df <- as.data.frame(pred_obj)
  names(pred_df)[1] <- "exposure_val"
  
  ref_x <- quantile(df[[exposure]], probs=ref_pct, na.rm=TRUE)
  
  if(family=="binomial") {
    ref_y <- pred_df$yhat[which.min(abs(pred_df$exposure_val-ref_x))]
    pred_df$yhat_rel <- pred_df$yhat / ref_y
    pred_df$lower_rel <- pred_df$lower / ref_y
    pred_df$upper_rel <- pred_df$upper / ref_y
    y_col <- "yhat_rel"; low_col <- "lower_rel"; up_col <- "upper_rel"
  } else {
    pred_df$yhat_rel <- pred_df$yhat
    pred_df$lower_rel <- pred_df$lower
    pred_df$upper_rel <- pred_df$upper
    y_col <- "yhat_rel"; low_col <- "lower_rel"; up_col <- "upper_rel"
  }
  
  p <- ggplot(pred_df, aes(x=exposure_val, y=.data[[y_col]])) +
    theme_minimal(base_size=13) +
    geom_ribbon(aes(ymin=.data[[low_col]], ymax=.data[[up_col]]), alpha=0.2) +
    geom_line(size=1.1) +
    geom_vline(xintercept=ref_x, linetype="dashed", color="grey40")
  
  if(family=="binomial") p <- p + geom_hline(yintercept=1, linetype="dotted", color="#A6DDCE")
  
  list(fit=fit, pred_df=pred_df, plot=p, ref_x=ref_x)
}


## -------------------------
## Step 5 - Higher-level wrapper: run_rcs_pipeline. Looping over exposures and saving outputs
## -------------------------
# data: data.frame
# outcome: name of outcome
# exposures: character vector of exposures to analyze
# covariates: covariates to adjust for
# family: "binomial" or "gaussian"
# nk: either a single integer (knots) or numeric vector of knot locations (for rms)
# out_prefix: file prefix
run_rcs_pipeline <- function(data, outcome, exposures, covariates = NULL,
                             family = c("binomial", "gaussian"), nk = 4,
                             ref_pct = 0.25, out_prefix = "rcs_out", nsim = 200) {
  family <- match.arg(family)
  results <- list()
  for (ex in exposures) {
    message("Running RCS for exposure: ", ex)
    res <- tryCatch({
      run_rcs_one(data = data, outcome = outcome, exposure = ex,
                  covariates = covariates, family = family, nk = nk,
                  ref_pct = ref_pct, nsim = nsim, out_prefix = out_prefix)
    }, error = function(e) {
      message("Failed for exposure ", ex, ": ", e$message)
      NULL
    })
    results[[ex]] <- res
  }
  return(results)
}

## -------------------------
## Step 6 - Example runs on simulated data
## -------------------------
# Logistic (binary) outcome with exposure env_A
run_rcs_one(
  data = example_df,
  outcome = "binary_outcome",
  exposure = "env_A",
  covariates = c("age_yrs", "bmi_kgm2", "male_flag"),
  family = "binomial",
  nk = 4,
  ref_pct = 0.25,
  nsim = 300
)
print(rcs_res_bin$plot)

# Continuous outcome example (env_A)
rcs_res_cont <- run_rcs_one(
  data = example_df,
  outcome = "continuous_outcome",
  exposure = "env_A",
  covariates = c("age_yrs", "bmi_kgm2", "male_flag"),
  family = "gaussian",
  nk = 4,
  ref_pct = 0.25,
  nsim = 300
)

rcs_res_cont$plot
