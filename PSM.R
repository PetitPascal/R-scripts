#-------------------------------------------------------------------------------
## Reproducible & generalisable propensity score matching (PSM) analysis pipeline
# - Uses MatchIt + cobalt + tableone for matching & diagnostics
# - Example uses 'lalonde' dataset (available in MatchIt)
# - Replace df_example and function arguments to run on your own data
#-------------------------------------------------------------------------------

## -------------------------
## Step 1 - Setup
## -------------------------
required <- c("MatchIt", "cobalt", "tableone", "dplyr", "ggplot2", "gridExtra", "openxlsx")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

set.seed(2025)

## -------------------------
## Step 2 - Using example dataset Lalonde available in MatchIt
## -------------------------

data("lalonde", package = "MatchIt")

df_example <- lalonde %>% 
  dplyr::rename(age_yrs = age,
    educ_years = educ,
    married_flag = married,
    no_degree = nodegree,
    inc_1974 = re74,
    inc_1975 = re75,
    inc_1978 = re78,
    treat_flag=treat)


## -------------------------
## Step 3 - Creating a helper function: check data, calculate SMDs and produce balance tables/plots
## -------------------------
psm_preflight <- function(df, treat_col, covs) {
  if (!treat_col %in% names(df)) stop("treatment column not found.")
  if (any(is.na(df[, treat_col]))) stop("NA in treatment column.")
  if (!all(df[[treat_col]] %in% c(0, 1))) {
    warning("Treatment column is not 0/1. Attempting to coerce (levels -> 0/1).")
    df[[treat_col]] <- as.numeric(as.factor(df[[treat_col]])) - 1
  }
  # reporting missingness
  miss <- colSums(is.na(df[, unique(c(treat_col, covs)), drop = FALSE]))
  message("Missingness per variable (0 = none):")
  print(miss[miss > 0 | TRUE][1:length(miss)]) # show all
  invisible(df)
}

smd_table <- function(mobj, covs, treat = "treat_flag", threshold = 0.1) {
  # Using cobalt::bal.tab to make SMD table
  btab <- cobalt::bal.tab(mobj, var.name = covs, un = TRUE, m.threshold = threshold)
  # Extracting numeric summary
  sumtab <- data.frame(
    variable = rownames(btab$Balance),
    SMD_un = btab$Balance[,"Diff.Un"],
    SMD_after = btab$Balance[,"Diff.Adj"]
  )
  sumtab <- sumtab[order(abs(sumtab$SMD_un), decreasing = TRUE), ]
  return(list(bal.tab = btab, smd = sumtab))
}

## -------------------------
## Step 4 - Main runner: run_psm_pipeline()
## -------------------------
# df: data.frame
# treat: name of treatment column (0/1)
# covariates: character vector of covariate names used to estimate PS
# method: "nearest", "optimal", "full", "genetic"
# ratio: matching ratio (1 = 1:1)
# caliper: caliper on the logit of PS (set NULL to skip)
# caliper_sensitivity: a numeric vector of calipers to test (optional)
# seed: RNG seed for reproducibility

run_psm_pipeline <- function(df, treat, covariates,
                             method = "nearest", ratio = 1, caliper = 0.2,
                             caliper_sensitivity = c(0.2, 0.1, 0.05),
                             distance = "logit", replace = FALSE,
                             std_caliper = TRUE, seed = 2025,
                             save_prefix = "psm_output") {
  set.seed(seed)
  psm_report <- list()
  
  # Step 1 - Preflight checks
  psm_preflight(df, treat, covariates)
  
  # Step 2 - Creating propensity score formula
  fmla <- as.formula(paste(treat, "~", paste(covariates, collapse = " + ")))
  message("Propensity score model formula: ", deparse(fmla))
  
  # Step 3 - Fitting an initial PS model (glm) and inspecting overlap
  ps_mod <- glm(fmla, data = df, family = binomial)
  df$pscore <- predict(ps_mod, type = "response")
  psm_report$ps_model <- ps_mod
  
  # Step 4 - Making an overlap plot (histograms by group)
  p_overlap <- ggplot2::ggplot(df, aes(x = pscore, fill = factor(.data[[treat]]))) +
    geom_histogram(position = "identity", alpha = 0.5, bins = 40) +
    labs(fill = "Treated", title = "Propensity score overlap", x = "PS", y = "Count") +
    theme_minimal()
  print(p_overlap)
  ggsave(paste0(save_prefix, "_ps_overlap.png"), p_overlap, width = 7, height = 4, dpi = 300)
  psm_report$ps_overlap_plot <- p_overlap
  
  # Step 5 - Running matchit (with or without caliper)
  match_args <- list(formula = fmla, data = df, method = method,
                     distance = distance, replace = replace, ratio = ratio)
  if (!is.null(caliper)) {
    match_args$caliper <- caliper
    match_args$std.caliper <- std_caliper
  }
  message("Running matchit with method='", method, "' and caliper=", ifelse(is.null(caliper), "NULL", caliper))
  m.out <- do.call(MatchIt::matchit, match_args)
  psm_report$matchit_obj <- m.out
  
  # Step 6 - Balance diagnostics (cobalt)
  bal <- cobalt::bal.tab(m.out, un = TRUE, m.threshold = 0.1)
  print(bal)  # console summary
  psm_report$bal_tab <- bal
  
  # Step 7 - Creating a love plot
  lp <- cobalt::love.plot(bal, stats = c("m"), abs = TRUE, var.order = "unadjusted",
                          threshold = 0.1, title = "Love plot: Absolute SMDs before & after")
  print(lp)
  
  # Step 8 - Saving the love plot
  ggsave(filename = paste0(save_prefix, "_loveplot.png"), plot = lp, width = 7, height = 6, dpi = 300)
  psm_report$love_plot <- lp
  
  # Step 9 - Extracting matched data and re-checking balance using tableone
  df_matched <- MatchIt::match.data(m.out)
  
  # Step 10 - Saving matched dataset
  write.csv(df_matched, paste0(save_prefix, "_matched_data.csv"), row.names = FALSE)
  message("Saved matched dataset to: ", paste0(save_prefix, "_matched_data.csv"))
  psm_report$matched_data <- df_matched
  
  # TableOne summary before/after
  tabVars <- covariates
  factorVars <- tabVars[sapply(df[, tabVars, drop = FALSE], function(x) is.factor(x) || is.character(x))]
  table_before <- tableone::CreateTableOne(vars = tabVars, strata = treat, data = df, factorVars = factorVars)
  table_after <- tableone::CreateTableOne(vars = tabVars, strata = treat, data = df_matched, factorVars = factorVars)
  psm_report$tableone_before <- table_before
  psm_report$tableone_after <- table_after
  print(table_before, smd = TRUE)
  print(table_after, smd = TRUE)
  
  # Step 11 - SMD table (detailed)
  st <- smd_table(m.out, covariates, treat = treat)
  psm_report$smd_table <- st$smd
  
  write.csv(st$smd, paste0(save_prefix, "_smd_table.csv"), row.names = FALSE) # writing SMD table to disk
  message("Saved SMD table to: ", paste0(save_prefix, "_smd_table.csv"))
  
  # Step 12 - Optional: caliper sensitivity analysis (if caliper_sensitivity provided)
  cal_results <- NULL
  if (!is.null(caliper_sensitivity) && length(caliper_sensitivity) > 0) {
    message("Running caliper sensitivity analysis for calipers: ", paste(caliper_sensitivity, collapse = ", "))
    cal_results <- data.frame(caliper = numeric(0), matched_n = integer(0), mean_smd_after = numeric(0))
    for (cval in caliper_sensitivity) {
      set.seed(seed)
      try({
        mtmp <- MatchIt::matchit(formula = fmla, data = df, method = method,
                                 distance = distance, caliper = cval, std.caliper = std_caliper, ratio = ratio, replace = replace)
        mm <- MatchIt::match.data(mtmp)
        baltmp <- cobalt::bal.tab(mtmp, un = TRUE)
        mean_smd <- mean(abs(baltmp$Balance[,"Diff.Adj"]), na.rm = TRUE)
        cal_results <- rbind(cal_results, data.frame(caliper = cval, matched_n = nrow(mm), mean_smd_after = mean_smd))
      }, silent = TRUE)
    }
    psm_report$caliper_sensitivity <- cal_results
    message("Caliper sensitivity results:\n")
    print(cal_results)
    write.csv(cal_results, paste0(save_prefix, "_caliper_sensitivity.csv"), row.names = FALSE)
    message("Saved caliper sensitivity to: ", paste0(save_prefix, "_caliper_sensitivity.csv"))
  }
  
  # Step 9 - Optional: estimating treatment effect on outcome (simple example)
  # Logistic regression of outcome ~ treat on matched sample (if outcome exists)
  # You can replace "outcome_name" with your outcome variable
  # Example: if 're78' exists in lalonde it will be used as outcome
  psm_report$notes <- "Use df_matched to run outcome models (e.g., glm on matched sample or using robust SEs)."
  
  return(psm_report)
}

## -------------------------
## Step 5 - Running example
## -------------------------
# Choose covariates (avoid re-using the original script's variable names)
covariate_list <- c("age_yrs", "educ_years", "race", "married_flag",
                    "no_degree", "inc_1974", "inc_1975", "inc_1978")

psm_out_example <- run_psm_pipeline(
  df = df_example,
  treat = "treat_flag",
  covariates = covariate_list,
  method = "nearest",
  ratio = 1,
  caliper = 0.2,
  caliper_sensitivity = c(0.25, 0.2, 0.1, 0.05),
  distance = "logit",
  save_prefix = "lalonde_psm"
)

## -------------------------
## Step 6 - How to use on your own data
## -------------------------
# Replace df_example with your own data frame (for example readxl::read_excel or openxlsx::read.xlsx)
# Example:
# df_user <- openxlsx::read.xlsx("path/your_file.xlsx", sheet = 1)
# psm_user_res <- run_psm_pipeline(df = df_user, treat = "exposure_col", covariates = c("age", "sex", "smoking"),
#                                  method = "nearest", ratio = 1, caliper = 0.2, save_prefix = "user_psm")
#
# After running, matched data is saved as user_psm_matched_data.csv and you can run outcome models on it:
# df_user_matched <- read.csv("user_psm_matched_data.csv")
# glm_outcome <- glm(outcome ~ treat_flag + covariates, data = df_user_matched, family = binomial)
