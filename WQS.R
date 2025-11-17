################################################################################
## Generalisable Weighted Quantile Sum regression (WQS) pipeline (self-contained)
## - Simulates correlated chemical mixture + covariates + binary outcome
## - Runs gWQS with sensible defaults and improvements
## - Provides helper function `run_wqs_pipeline()` for user data
## - Saves outputs: weights CSV, model object, plots
################################################################################

## -------------------------
## Step 1 - Setup
## -------------------------
rm(list = ls())
required_pkgs <- c("gWQS", "openxlsx", "MASS", "dplyr", "ggplot2", "RColorBrewer", "tidyr", "forcats")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
}

set.seed(2025)   # reproducibility seed (change as needed)

## -------------------------
## Step 2 - Creating utility helper functions
## -------------------------

# Safe log1p transform that handles negative/NA
safe_log1p <- function(x) {
  x <- ifelse(is.na(x), NA_real_, x)
  x[x < 0] <- 0  # chemical concentrations should be >=0; if negative due to data error, clip to 0
  log1p(x)
}

# Function to plot WQS weights with threshold line
plot_wqs_weights <- function(weights_df, top_n = NULL, title = "WQS: final weights") {
  df <- weights_df
  if (!is.null(top_n)) df <- df %>% arrange(desc(mean_weight)) %>% head(top_n)
  df <- df %>% mutate(mix = forcats::fct_reorder(mix, mean_weight))
  p <- ggplot(df, aes(x = mix, y = mean_weight)) +
    geom_col(fill = "#A6DDCE", color = "black") +
    coord_flip() +
    geom_hline(yintercept = 1 / nrow(df), linetype = "dashed", color = "red") +
    labs(x = "", y = "Mean weight", title = title, subtitle = paste0("Dashed line = 1 / K (K = ", nrow(df), ")")) +
    theme_minimal()
  return(p)
}

## -------------------------
## Step 3 - Simulating a dataset (replace or skip this block to use your own data)
## -------------------------
n <- 800
K <- 14  # number of mixture components

# Creating correlated exposures via MVN then exponentiate to look like concentrations
rho <- 0.55
Sigma <- rho ^ as.matrix(dist(1:K))
raw <- MASS::mvrnorm(n = n, mu = rep(0, K), Sigma = Sigma)

# Making positive concentration-like variables (scaling down to avoid huge range)
exposures <- as.data.frame(exp(raw / 3))
exposure_names <- paste0("mix_", seq_len(K))
colnames(exposures) <- exposure_names

# Covariates
covariates <- data.frame(
  age_years = round(rnorm(n, mean = 50, sd = 12)),
  education_high = rbinom(n, 1, prob = 0.35), # binary indicator for "higher education"
  residence_urban = rbinom(n, 1, prob = 0.6)
)

# Simulating true mixture effect (non-linear + interactions)
true_weights <- c(0.4, 0.2, rep(0.05, K - 2))
h_z <- as.numeric(as.matrix(log1p(exposures)) %*% true_weights)
h_z <- h_z + 0.6 * (log1p(exposures$mix_1) * log1p(exposures$mix_2))

# Binary outcome probability
linpred <- -2 + 0.02 * covariates$age_years + 0.35 * covariates$education_high + 0.3 * covariates$residence_urban + h_z
prob_case <- 1 / (1 + exp(-linpred))
case_var <- rbinom(n, 1, prob_case)

# Combining into one data.frame
sim_data <- dplyr::bind_cols(
  tibble::tibble(subject_id = seq_len(n), case_status = case_var),
  covariates,
  exposures
)

# Quick check
message("Simulated data: ", nrow(sim_data), " rows; ", length(exposure_names), " exposures.")

## -------------------------
## Step 4 - Preprocessing suggestions & quick EDA
## -------------------------

# Checking missingness
na_summary <- colSums(is.na(sim_data))
print("NA counts (first 20 variables):")
print(na_summary[1:min(20, length(na_summary))])

# Distribution checks for exposures (optional)
boxplot(log1p(sim_data[, exposure_names]), main = "log1p exposures (simulated)")

## -------------------------
## Step 5 - Running a single WQS model
## -------------------------
# Important gWQS notes:
# - mix_name: vector of names of mixture variables in data
# - formula: outcome ~ wqs + covariates
# - q: number of quantiles (3 or 4 typical)
# - validation: proportion used for validation (0.4-0.7 common)
# - b: number bootstrap samples for weight estimation (500-2000 common)
# - b1_pos: TRUE if we expect positive association (otherwise FALSE)
# - family: "binomial" for binary outcome; "gaussian" for continuous

# Choosing settings
q_val <- 4
validation_frac <- 0.6
n_boot <- 1000
b1_positive <- TRUE
direction_to_test <- if (b1_positive) "positive" else "negative"

message("Running gWQS with q=", q_val, ", validation=", validation_frac, ", b=", n_boot, ", direction=", direction_to_test)

# A small wrapper to run gWQS safely (with log transform option)
run_single_gwqs <- function(data, outcome_name, mix_vars, covariate_names = NULL,
                            q = 4, validation = 0.6, b = 1000, b1_pos = TRUE,
                            family = "binomial", seed = 2025, log_transform = TRUE,
                            table_output_path = "wqs_weights.csv") {
  
  set.seed(seed)
  df <- data
  
  # ---- Optional log transform ----
  if (log_transform) {
    df[, mix_vars] <- lapply(df[, mix_vars, drop = FALSE], function(x) {
      x[x < 0] <- 0      # avoid negative chemical values
      return(log1p(x))   # safe log transform
    })
  }
  
  # ---- Building model formula ----
  cov_part <- if (!is.null(covariate_names) && length(covariate_names) > 0) {
    paste(covariate_names, collapse = " + ")
  } else {
    ""
  }
  
  model_formula <- if (cov_part == "") {
    as.formula(paste(outcome_name, "~ wqs"))
  } else {
    as.formula(paste(outcome_name, "~ wqs +", cov_part))
  }
  
  # ---- Binary-family check ----
  if (family == "binomial") {
    if (!all(df[[outcome_name]] %in% c(0, 1))) {
      stop("For family='binomial', the outcome must be coded 0/1.")
    }
  }
  
  # ---- Running gWQS ----
  gwqs_out <- gWQS::gwqs(
    formula = model_formula,
    mix_name = mix_vars,
    data = df,
    q = q,
    validation = validation,
    b = b,
    b1_pos = b1_pos,
    family = family,
    seed = seed,
    plots = FALSE,
    tables = FALSE
  )
  
  # ---- Extracting final weights ----
  w <- gwqs_out$final_weights
  
  # Automatic renaming depending on number of columns
  if (ncol(w) == 2) {
    colnames(w) <- c("mix", "mean_weight")
  } else if (ncol(w) == 3) {
    colnames(w) <- c("mix", "mean_weight", "sd_weight")
  } else if (ncol(w) == 4) {
    colnames(w) <- c("mix", "mean_weight", "sd_weight", "median_weight")
  } else {
    stop("Unexpected number of columns in gwqs_out$final_weights: ", ncol(w))
  }
  
  w <- w[order(-w$mean_weight), ]
  
  write.csv(w, file = table_output_path, row.names = FALSE)
  message("Saved final weights to: ", normalizePath(table_output_path))
  
  return(list(gwqs = gwqs_out, final_weights = w))
}

# Running the WQS on simulated data
wqs_res <- run_single_gwqs(
  data = sim_data,
  outcome_name = "case_status",
  mix_vars = exposure_names,
  covariate_names = c("age_years", "education_high", "residence_urban"),
  q = q_val,
  validation = validation_frac,
  b = n_boot,
  b1_pos = b1_positive,
  family = "binomial",
  seed = 2025,
  log_transform = TRUE,
  table_output_path = "wqs_sim_final_weights.csv"
)

# Viewing results
print(summary(wqs_res$gwqs))
print("Top weights:")
print(head(wqs_res$final_weights, 10))

# Plotting weights (full)
p_w <- plot_wqs_weights(wqs_res$final_weights, title = "Simulated data: WQS final weights")+
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
print(p_w)
ggsave("wqs_weights_full.png", p_w, width = 8, height = 6, dpi = 300)

## -------------------------
## Step 6 - Additional visualization: barplot for top-K contributors
## -------------------------
topK <- 12
p_topk <- plot_wqs_weights(wqs_res$final_weights, top_n = topK, title = paste0("Top ", topK, " mixture weights"))+
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
print(p_topk)
ggsave("wqs_weights_topK.png", p_topk, width = 8, height = 6, dpi = 300)

## -------------------------
## Step 7 - Sensitivity analysis — repeated holdout WQS (simple manual loop)
## Purpose: evaluate stability of estimated weights across repeated random splits
## -------------------------
run_repeated_wqs <- function(data, outcome_name, mix_vars, covariate_names = NULL,
                             q = 4, validation = 0.6, b = 500, b1_pos = TRUE,
                             family = "binomial", reps = 10, log_transform = TRUE, seed = 123) {
  set.seed(seed)
  all_weights <- vector("list", reps)
  for (i in seq_len(reps)) {
    seed_i <- seed + i
    res_i <- run_single_gwqs(
      data = data,
      outcome_name = outcome_name,
      mix_vars = mix_vars,
      covariate_names = covariate_names,
      q = q,
      validation = validation,
      b = b,
      b1_pos = b1_pos,
      family = family,
      seed = seed_i,
      log_transform = log_transform,
      table_output_path = tempfile(fileext = ".csv")
    )
    all_weights[[i]] <- res_i$final_weights %>% mutate(rep = i)
    message("Completed repeated WQS run ", i, "/", reps)
  }
  all_w_df <- do.call(rbind, all_weights)
  return(all_w_df)
}

# Running modest repeated holdout (can increase reps but will take longer)
reps <- 6
message("Running repeated-holdout sensitivity with ", reps, " repeats (b=", 400, ") ...")
repeated_w <- run_repeated_wqs(
  data = sim_data,
  outcome_name = "case_status",
  mix_vars = exposure_names,
  covariate_names = c("age_years", "education_high", "residence_urban"),
  q = q_val,
  validation = validation_frac,
  b = 400,
  b1_pos = b1_positive,
  family = "binomial",
  reps = reps,
  log_transform = TRUE,
  seed = 2025
)

# Summarizing variation in weights across repeats
weight_summary <- repeated_w %>%
  group_by(mix) %>%
  summarize(mean_weight = mean(mean_weight), sd_weight = sd(mean_weight), n_reps = n()) %>%
  arrange(desc(mean_weight))

print("Summary of repeated-holdout weights (first 10):")
print(head(weight_summary, 10))

# Plotting mean +/- sd across repeats
p_rep <- ggplot(weight_summary %>% head(20), aes(x = forcats::fct_reorder(mix, mean_weight), y = mean_weight)) +
  geom_col(fill = "#8DD3C7", color = "black") +
  geom_errorbar(aes(ymin = mean_weight - sd_weight, ymax = mean_weight + sd_weight), width = 0.3) +
  coord_flip() +
  labs(x = "", y = "Mean weight (across repeats)", title = "WQS repeated-holdout: mean ± sd of weights") +
  theme_minimal()
print(p_rep)
ggsave("wqs_repeated_holdout_weights.png", p_rep, width = 9, height = 6, dpi = 300)

## -------------------------
## Step 8 - Saving full model object for reproduction
## -------------------------
saveRDS(wqs_res$gwqs, file = "wqs_model_sim.rds")
message("Saved WQS model object to wqs_model_sim.rds")

## -------------------------
## Creating the final helper function
## -------------------------
# Usage: supply your dataframe, outcome name, mixture variable names, covariate names

# Example
user_res <- run_wqs_pipeline(user_df, outcome_col = "y", mix_cols = my_chem_names,
                             covariates = c("age","sex"), family = "binomial",
                             q = 4, validation = 0.6, b = 1000, reps = 20)

run_wqs_pipeline <- function(data, outcome_col, mix_cols, covariates = NULL,
                             q = 4, validation = 0.6, b = 1000, b1_pos = TRUE,
                             family = "binomial", log_transform = TRUE, reps = 0,
                             out_prefix = "wqs_user", seed = 2025) {
  # Basic checks
  if (!all(mix_cols %in% names(data))) stop("Some mixture columns are missing from data.")
  if (!outcome_col %in% names(data)) stop("Outcome column not found.")
  if (family == "binomial" && !all(data[[outcome_col]] %in% c(0, 1))) stop("Binary outcome must be 0/1.")
  
  # Running a single model and saving
  single_res <- run_single_gwqs(
    data = data, outcome_name = outcome_col, mix_vars = mix_cols,
    covariate_names = covariates, q = q, validation = validation, b = b,
    b1_pos = b1_pos, family = family, seed = seed, log_transform = log_transform,
    table_output_path = paste0(out_prefix, "_final_weights.csv")
  )
  saveRDS(single_res$gwqs, file = paste0(out_prefix, "_model.rds"))
  
  # Optionally running repeated holdouts for sensitivity
  if (reps > 0) {
    repeated_df <- run_repeated_wqs(
      data = data, outcome_name = outcome_col, mix_vars = mix_cols, covariate_names = covariates,
      q = q, validation = validation, b = round(b / 2), b1_pos = b1_pos, family = family, reps = reps,
      log_transform = log_transform, seed = seed
    )
    write.csv(repeated_df, file = paste0(out_prefix, "_repeated_weights.csv"), row.names = FALSE)
    message("Saved repeated-holdout weights to: ", paste0(out_prefix, "_repeated_weights.csv"))
  } else {
    repeated_df <- NULL
  }
  
  return(list(single = single_res, repeated = repeated_df))
}

## -------------------------
## Example of use on simulated data
## -------------------------

user_example <- run_wqs_pipeline(sim_data, outcome_col = "case_status",
                                 mix_cols = exposure_names,
                                 covariates = c("age_years","education_high","residence_urban"),
                                 q = 4, validation = 0.6, b = 1000, reps = 10,
                                 out_prefix = "sim_wqs", seed = 2025)

cat("\nWQS pipeline complete. Files produced:\n")
cat("- wqs_sim_final_weights.csv (weights table)\n")
cat("- wqs_weights_full.png, wqs_weights_topK.png (plots)\n")
cat("- wqs_repeated_holdout_weights.png (if repeated run was used)\n")
cat("- wqs_model_sim.rds (model object)\n")
cat("\nUse run_wqs_pipeline() to run on your own dataset — see comments above.\n")
