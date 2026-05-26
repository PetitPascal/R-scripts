#-------------------------------------------------------------------------------
## Reproducible & Generalisable Poisson / Negative Binomial Regression Script
# Covers:
#   - Poisson regression (simple and multiple)
#   - Overdispersion testing and correction (quasi Poisson, negative binomial)
#   - Zero-inflation checks
#   - Assumption testing: distribution fit, residuals, influential obs
#   - Rate models with offset
#   - Incidence Rate Ratio (IRR) extraction and plots
#   - Reusable pipeline function
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup
rm(list = ls())
set.seed(123)

required_pkgs <- c("MASS", "dplyr", "tidyr", "ggplot2", "broom","gtsummary", "car", "lmtest", "performance",
                   "pscl", "AER", "ggeffects", "forestmodel","ggstats", "tibble", "conflicted")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# for avoiding package conflicts
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("tidy","broom")
conflicted::conflicts_prefer(stats::predict)

#--------------------------------------------
## Step 2: Simulate datasets
n<-500

# Predictors
age<-round(rnorm(n, 50, 12))
sex<-factor(rbinom(n, 1, 0.48), labels = c("Female","Male"))
bmi<-rnorm(n, 26, 4)
smoking<-factor(rbinom(n, 1, 0.3),  labels = c("No","Yes"))
exposure<-rnorm(n, 5, 2)
country<-factor(sample(c("FR","DE","IT","ES"), n, replace = TRUE))

# Log-linear predictor
log_mu<-1.5 + 0.02*age + 0.3*(sex=="Male") +0.04*bmi + 0.5*(smoking=="Yes") + 0.15*exposure

# Dataset 1: Poisson (no overdispersion)
y_pois<-rpois(n, exp(log_mu))

# Dataset 2: Overdispersed (negative binomial, theta = 2)
y_nb<-MASS::rnegbin(n, mu = exp(log_mu), theta = 2)

# Dataset 3: With population offset (rate model)
log_pop<-log(round(runif(n, 1000, 100000)))
y_rate<-rpois(n, exp(log_mu + log_pop - mean(log_mu + log_pop) + 3))

sim_data<-tibble(y_pois= y_pois,
  y_nb=y_nb,
  y_rate=y_rate,
  log_pop=log_pop,
  age=age,
  sex=sex,
  bmi=bmi,
  smoking=smoking,
  exposure=exposure,
  country=country)

#--------------------------------------------
## Step 3: EDA

# Distribution of count outcomes
par(mfrow = c(1,2))
hist(y_pois, main = "Poisson outcome", col = "#A6DDCE", xlab = "Count")
hist(y_nb,   main = "Neg. Binomial outcome", col = "#d6604d", xlab = "Count")
par(mfrow = c(1,1))

# Theoretical vs empirical Poisson distribution check
theo_pois<-rpois(n, mean(y_pois))
dist_check<-tibble(value= c(y_pois, theo_pois),
  distribution = rep(c("Empirical","Theoretical"), each = n))

ggplot(dist_check, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Poisson: empirical vs theoretical distribution",
       x = "Count", y = "Density") +
  scale_fill_manual(values = c("#2166ac","#d6604d")) +
  theme_bw(base_size = 14)

#--------------------------------------------
## Simple Poisson regression

mod_pois_simple<-glm(y_pois ~ exposure,
                       data=sim_data,
                       family=poisson(link = "log"))
summary(mod_pois_simple)

# Incidence Rate Ratio
irr_simple<-broom::tidy(mod_pois_simple,conf.int = TRUE, exponentiate = TRUE)
irr_simple

# Plot
ggplot(sim_data, aes(x = exposure, y = y_pois)) +
  geom_point(alpha = 0.3) +
  geom_line(aes(y = fitted(mod_pois_simple)), color = "#2166ac",
            linewidth = 1.2) +
  labs(title = "Simple Poisson: observed (points) vs fitted (line)",
       x = "Exposure", y = "Count") +
  theme_bw(base_size = 14)

#--------------------------------------------
## Multiple Poisson regression

# Fitting model
mod_pois<-glm(y_pois ~ exposure + age + sex + bmi + smoking,
                data   = sim_data,
                family = poisson(link = "log"))
summary(mod_pois)

# Assumption checks
check_poisson_assumptions<-function(model){
  
  # Goodness of fit (chi-squared)
  gof_stat<-model$deviance
  gof_df<-model$df.residual
  gof_p<-pchisq(gof_stat, gof_df, lower.tail = FALSE)
  cat("Goodness-of-fit chi2:", round(gof_stat, 3),
      "| df:", gof_df, "| p:", round(gof_p, 4),
      "(p > 0.05 = adequate fit)\n")
  
  # Overdispersion (deviance / df)
  disp<-gof_stat/gof_df
  cat("Overdispersion (dev/df):", round(disp, 3),
      if (disp > 1.5) "*** OVERDISPERSION — use quasiPoisson or NB ***"
      else "OK", "\n")
  
  # Formal overdispersion test
  if (requireNamespace("AER", quietly = TRUE)) {
    dt <- AER::dispersiontest(model)
    cat("Dispersion test: z =", round(dt$statistic, 3),
        "| p =", round(dt$p.value, 4),
        "(p < 0.05 = overdispersion)\n")
  }
  
  # Performance package check
  perf_od <- performance::check_overdispersion(model)
  print(perf_od)
  
  # Zero inflation
  obs_zeros <- sum(model$y == 0)
  exp_zeros <- sum(dpois(0, fitted(model)))
  cat("Observed zeros:", obs_zeros,
      "| Expected zeros:", round(exp_zeros, 1),
      if (obs_zeros > 2 * exp_zeros) "*** ZERO INFLATION ***" else "OK",
      "\n")
  
  # Multicollinearity
  vif_r <- tryCatch(car::vif(model), error = function(e) NULL)
  if (!is.null(vif_r))
    cat("Max VIF:", round(max(vif_r), 3),
        if (max(vif_r) > 5) "*** HIGH VIF ***" else "OK", "\n")
  
  # Influential observations
  cooks  <- cooks.distance(model)
  n_infl <- sum(cooks > 4 / length(model$y))
  cat("Influential obs (Cook > 4/n):", n_infl, "\n")
  
  # Residual plots
  par(mfrow = c(2,2))
  plot(model)
  par(mfrow = c(1,1))
  
  invisible(list(disp = disp, gof_p = gof_p,
                 zeros = obs_zeros, exp_zeros = exp_zeros))
}

pois_checks <- check_poisson_assumptions(mod_pois)

# IRR table and plot
extract_irr<-function(model){
  broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(sig=case_when(
        conf.low  > 1 ~ "IRR > 1",
        conf.high < 1 ~ "IRR < 1",
        TRUE          ~ "NS"),
      sig = factor(sig, levels = c("IRR < 1","NS","IRR > 1")))
}

plot_irr<-function(irr_df, title = "Incidence Rate Ratios (95% CI)") {
  ggplot(irr_df, aes(x = reorder(term, estimate),
                     y = estimate,
                     ymin = conf.low, ymax = conf.high,
                     color = sig)) +
    geom_pointrange(size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "darkblue") +
    coord_flip() +
    scale_color_manual(
      values = c("IRR < 1" = "darkgreen",
                 "NS"      = "black",
                 "IRR > 1" = "red"),
      name = "") +
    labs(title = title, x = "", y = "IRR (95% CI)") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")
}

irr_pois <- extract_irr(mod_pois)
irr_pois
plot_irr(irr_pois, "Poisson: Incidence Rate Ratios")

# Summary table
mod_pois %>%
  gtsummary::tbl_regression(exponentiate = TRUE) %>%
  gtsummary::bold_p()

# Forest plot
forestmodel::forest_model(mod_pois, exponentiate = TRUE)

# Predicted counts plot
plot(ggeffects::ggeffect(mod_pois, terms = "exposure")) +
  labs(title = "Poisson: predicted count by exposure")

plot(ggeffects::ggeffect(mod_pois, terms = c("exposure","smoking"))) +
  labs(title = "Predicted count: exposure × smoking")

#--------------------------------------------
## Quasi Poisson regression - use when deviance/df > 1.5 (overdispersion present)

mod_quasi<-glm(y_nb ~ exposure + age + sex + bmi + smoking,
                 data= sim_data,
                 family= quasipoisson(link = "log"))
summary(mod_quasi)

#  Type II effects
car::Anova(mod_quasi, test = "F")

# IRR
broom::tidy(mod_quasi, conf.int = TRUE, exponentiate = TRUE)

#--------------------------------------------
## Negative binomial regression: better than quasi Poisson when variance >> mean (count data)

mod_nb<-MASS::glm.nb(y_nb ~ exposure + age + sex + bmi + smoking,data = sim_data)
summary(mod_nb)

cat("\nTheta (NB dispersion):", round(mod_nb$theta, 4),"(lower = more overdispersion)\n")

# IRR
irr_nb<-extract_irr(mod_nb)
irr_nb
plot_irr(irr_nb, "Negative Binomial: IRR")

# LRT: Poisson vs. NB (test for overdispersion)
lrt_nb<-lmtest::lrtest(mod_pois, mod_nb)
lrt_nb

#--------------------------------------------
## Rate model with offset: use when observations represent different population sizes/time periods
# offset = log(population) to model rates rather than raw counts

mod_rate<-glm(y_rate ~ exposure + age + sex + smoking + offset(log_pop),
                data=sim_data,
                family=poisson(link = "log"))
summary(mod_rate)

# IRR
broom::tidy(mod_rate, conf.int = TRUE, exponentiate = TRUE)

# Overdispersion check for rate model
disp_rate<-mod_rate$deviance / mod_rate$df.residual
disp_rate

if (disp_rate > 1.5) {
  cat("Overdispersion detected — fitting NB rate model.\n")
  mod_rate_nb<-MASS::glm.nb(y_rate ~ exposure + age + sex + smoking + offset(log_pop),data=sim_data)
  summary(mod_rate_nb)
}

#--------------------------------------------
## Zero inflated Poisson (ZIP): use when there are excess zeros beyond what Poisson predicts

# Simulating zero-inflated data
prob_zero<-plogis(-1 + 0.3*(sim_data$smoking=="Yes"))
zi_flag<-rbinom(n, 1, prob_zero)
y_zip<-ifelse(zi_flag == 1, 0, rpois(n, exp(log_mu)))
sim_data$y_zip<-y_zip

# Fitting ZIP model
mod_zip<-pscl::zeroinfl(y_zip ~ exposure + age + sex + bmi + smoking |   # count part
         smoking + age, # zero-inflation part
  data = sim_data, dist = "poisson")
summary(mod_zip)

# Zero-inflated NB
mod_zinb<-pscl::zeroinfl(y_zip ~ exposure + age + sex + bmi + smoking |
    smoking + age,
  data = sim_data, dist = "negbin")
summary(mod_zinb)

# Comparing ZIP vs. Zero-inflated NB via Vuong test
vuong_res<-pscl::vuong(mod_zip, mod_zinb)
vuong_res

# AIC comparison
aic_df <- data.frame(
  Model= c("Poisson","Quasi Poisson","NB","ZIP","ZINB"),
  AIC= c(AIC(mod_pois), NA, AIC(mod_nb), AIC(mod_zip), AIC(mod_zinb)))
aic_df

#--------------------------------------------
## Model comparison

mod_pois_null<-glm(y_pois ~ 1, data = sim_data, family = poisson())

# LRT
print(anova(mod_pois_null, mod_pois, test = "Chisq"))

# AIC/BIC
AIC(mod_pois_null, mod_pois, mod_nb) %>% arrange(AIC)

# Summary table
gtsummary::tbl_merge(
  tbls = list(gtsummary::tbl_regression(mod_pois, exponentiate = TRUE),
    gtsummary::tbl_regression(mod_nb,   exponentiate = TRUE)),
  tab_spanner = c("Poisson","Negative Binomial"))

# ggstats comparison
ggstats::ggcoef_compare(list("Poisson" = mod_pois, "Neg. Binomial" = mod_nb),
  type= "faceted",
  intercept= FALSE,
  exponentiate=TRUE) + 
  labs(title= "Poisson vs NB: coefficient comparison")

#--------------------------------------------
## Reusable pipeline
run_count_pipeline<-function(data, outcome_col,
                             predictor_cols,
                             offset_col=NULL,
                             family= c("poisson","nb","zip","zinb"),
                             zi_predictors=NULL,
                             seed=123) {
  set.seed(seed)
  family<-match.arg(family)
  df<-data
  
  pred_str<-paste(predictor_cols, collapse = " + ")
  offset_str<-if (!is.null(offset_col))
    paste0(" + offset(", offset_col, ")") else ""
  
  if(family %in% c("poisson","nb")) {
    form<-as.formula(paste(outcome_col, "~", pred_str, offset_str))
    fit<-if (family == "poisson")
      glm(form,  data = df, family = poisson())
    else
      MASS::glm.nb(form, data = df)
  } else {
    zi_str<-if (!is.null(zi_predictors))
      paste(zi_predictors, collapse = " + ") else "1"
    form<-as.formula(paste(outcome_col, "~", pred_str, "|", zi_str))
    dist<-if (family == "zip") "poisson" else "negbin"
    fit<-pscl::zeroinfl(form, data = df, dist = dist)
  }
  
  cat("\nCount model:", family, "\n")
  print(summary(fit))
  
  # Overdispersion
  if (family == "poisson") {
    disp <- fit$deviance / fit$df.residual
    cat("Overdispersion:", round(disp, 3),
        if (disp > 1.5) "— consider NB or quasiPoisson" else "— OK", "\n")
  }
  
  # IRR
  if (family %in% c("poisson","nb")) {
    irr_df <- broom::tidy(fit, conf.int=TRUE, exponentiate=TRUE) %>%
      filter(term != "(Intercept)")
    cat("\nIRR\n"); print(irr_df)
    print(plot_irr(irr_df %>%
                     mutate(sig = case_when(conf.low>1~"IRR > 1",
                                            conf.high<1~"IRR < 1",
                                            TRUE~"NS"),
                            sig = factor(sig, c("IRR < 1","NS","IRR > 1"))),
                   paste("IRR:", family)))
  }
  
  return(fit)
}

# Examples
fit_pois_pipe<-run_count_pipeline(data= sim_data,
  outcome_col= "y_pois",
  predictor_cols=c("exposure","age","sex","bmi","smoking"),
  family="poisson")

fit_nb_pipe<-run_count_pipeline(data=sim_data,
  outcome_col="y_nb",
  predictor_cols=c("exposure","age","sex","bmi","smoking"),
  family="nb")

fit_zip_pipe<-run_count_pipeline(
  data=sim_data,
  outcome_col= "y_zip",
  predictor_cols=c("exposure","age","sex","bmi","smoking"),
  family="zip",
  zi_predictors=c("smoking","age"))
