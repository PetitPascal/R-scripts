#-------------------------------------------------------------------------------
## Reproducible & Generalisable Linear Regression Script
# Covers: simple LR, multiple LR, robust LR, elastic net, LASSO, ridge, partial least squares (PLS)
# - Simulates data with correlated predictors + continuous outcome
# - Tests all key assumptions for each method
# - Compares predictive performance via cross-validation
# - Includes run_linear_pipeline() wrapper
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup
rm(list = ls())
set.seed(123)

required_pkgs <- c("MASS", "dplyr", "ggplot2", "tidyr", "broom","glmnet", "caret", "pls", "spls", "car",
                   "lmtest", "performance", "sandwich","gtsummary", "corrplot")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

#--------------------------------------------
## Step 2: Simulating data

n<-500
p<-15 # predictors

# Correlated predictors
rho<-0.5
Sigma<-rho ^ as.matrix(dist(1:p))
X_raw<-MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
colnames(X_raw) <- paste0("X", seq_len(p))

# True coefficients (sparse: only first 5 non-zero)
beta_true<-c(1.5, -1.2, 0.8, -0.6, 0.4, rep(0, p - 5))
y<-X_raw %*% beta_true + rnorm(n, 0, 2)

sim_data<-as.data.frame(cbind(y = as.numeric(y), X_raw))
X_mat<-as.matrix(sim_data[, -1])
y_vec<-sim_data$y

summary(y_vec)

#--------------------------------------------
## Step 3: Assumption helper function

check_lm_assumptions <- function(model, label = "") {
  cat("\n========== Assumption checks:", label, "==========\n")
  
  # Normality of residuals
  sw  <- shapiro.test(residuals(model))
  cat("Shapiro-Wilk normality: W =", round(sw$statistic, 4),
      "| p =", round(sw$p.value, 4))
  if (sw$p.value < 0.05) cat(" *** NON-NORMAL ***")
  cat("\n")
  
  # Homoscedasticity (Breusch-Pagan)
  bp  <- lmtest::bptest(model)
  cat("Breusch-Pagan homoscedasticity: BP =", round(bp$statistic, 4),
      "| p =", round(bp$p.value, 4))
  if (bp$p.value < 0.05) cat(" *** HETEROSCEDASTIC ***")
  cat("\n")
  
  # Independence (Durbin-Watson)
  dw  <- lmtest::dwtest(model)
  cat("Durbin-Watson autocorrelation: DW =", round(dw$statistic, 4),
      "| p =", round(dw$p.value, 4))
  if (dw$p.value < 0.05) cat(" *** AUTOCORRELATION ***")
  cat("\n")
  
  # Multicollinearity (VIF)
  vif_res <- tryCatch(car::vif(model), error = function(e) NULL)
  if (!is.null(vif_res)) {
    cat("Max VIF:", round(max(vif_res), 3))
    if (max(vif_res) > 5) cat(" *** HIGH MULTICOLLINEARITY (>5) ***")
    cat("\n")
  }
  
  # Influential observations (Cook's distance)
  cooks  <- cooks.distance(model)
  n_infl <- sum(cooks > 4 / length(cooks))
  cat("Influential obs (Cook's D > 4/n):", n_infl, "\n")
  
  # Diagnostic plots
  par(mfrow = c(2, 2))
  plot(model, main = label)
  par(mfrow = c(1, 1))
}

#--------------------------------------------
## Step 4: model fitting

#- - - - - -
# Simple linear regression
slr_fit <- lm(y ~ X1, data = sim_data)
summary(slr_fit)
broom::tidy(slr_fit)
check_lm_assumptions(slr_fit)

ggplot(sim_data, aes(x = X1, y = y)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, color = "#2166ac") +
  labs(title = "Simple linear regression: y ~ X1",
       x = "X1", y = "y") +
  theme_bw(base_size = 14)

#- - - - - -
# Multiple linear regression
mlr_fit <- lm(y ~ ., data = sim_data)
summary(mlr_fit)
check_lm_assumptions(mlr_fit, "Multiple LR")

tidy_mlr <- broom::tidy(mlr_fit, conf.int = TRUE)
tidy_mlr

ggplot(tidy_mlr %>% filter(term != "(Intercept)"),
       aes(x = reorder(term, estimate), y = estimate,
           ymin = conf.low, ymax = conf.high)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Multiple LR: coefficients (95% CI)",
       x = "", y = "Estimate") +
  theme_bw(base_size = 13)

#- - - - - -
# Robust linear regression
rlm_fit <- MASS::rlm(y ~ ., data = sim_data, method = "MM")
summary(rlm_fit)

# Using sandwich/lmtest for robust inference on rlm
rlm_coeftest <- lmtest::coeftest(rlm_fit,vcov. = sandwich::vcovHC(rlm_fit, type="HC3"))
rlm_coeftest

# R2
ss_res<-sum(residuals(rlm_fit)^2)
ss_tot<-sum((y_vec - mean(y_vec))^2)
r2_rlm<-1 - ss_res / ss_tot
cat("Pseudo-R2 (robust LR):", round(r2_rlm, 4), "\n")

# Comparing OLS vs. robust estimates for influential observations
sim_data$ols_resid<-residuals(mlr_fit)
sim_data$rlm_wts<-rlm_fit$w # Huber/bisquare weights (low = influential)

ggplot(sim_data, aes(x = seq_len(n), y = rlm_wts)) +
  geom_point(alpha = 0.5, color = "#d6604d") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(title = "Robust LR: observation weights (low = influential/outlier)",
       x = "Observation", y = "Weight") +
  theme_bw(base_size = 14)

#- - - - - -
# Ridge regression
cv_ridge <- glmnet::cv.glmnet(X_mat, y_vec, alpha = 0, nfolds = 10)
lambda_ridge_min <- cv_ridge$lambda.min
lambda_ridge_1se <- cv_ridge$lambda.1se
cat("Ridge: lambda.min =", round(lambda_ridge_min, 4),
    "| lambda.1se =", round(lambda_ridge_1se, 4), "\n")

plot(cv_ridge)

ridge_fit <- glmnet::glmnet(X_mat, y_vec, alpha = 0,lambda = lambda_ridge_1se)
ridge_coef <- as.matrix(coef(ridge_fit))
cat("\nRidge coefficients (lambda.1se)\n")
round(ridge_coef, 4)

# Coefficient path
plot(glmnet::glmnet(X_mat, y_vec, alpha = 0),
     xvar = "lambda", label = TRUE,
     main = "Ridge: coefficient paths")

#- - - - - -
# LASSO regression

cv_lasso <- glmnet::cv.glmnet(X_mat, y_vec, alpha = 1, nfolds = 10)
lambda_lasso_min <- cv_lasso$lambda.min
lambda_lasso_1se <- cv_lasso$lambda.1se
cat("LASSO: lambda.min =", round(lambda_lasso_min, 4),
    "| lambda.1se =", round(lambda_lasso_1se, 4), "\n")

plot(cv_lasso)

lasso_fit  <- glmnet::glmnet(X_mat, y_vec, alpha = 1,
                             lambda = lambda_lasso_1se)
lasso_coef <- as.matrix(coef(lasso_fit))
cat("\nLASSO coefficients (lambda.1se)\n")
round(lasso_coef, 4)

selected_vars_lasso <- rownames(lasso_coef)[lasso_coef != 0 &
                                              rownames(lasso_coef) != "(Intercept)"]
cat("LASSO-selected variables:", paste(selected_vars_lasso, collapse=", "), "\n")

# Coefficient path
plot(glmnet::glmnet(X_mat, y_vec, alpha = 1),
     xvar = "lambda", label = TRUE,
     main = "LASSO: coefficient paths")

#- - - - - -
# Elastic net regression

# Grid search over alpha
alpha_grid <- seq(0.1, 0.9, by = 0.1)
en_results <- lapply(alpha_grid, function(a) {
  cv <- glmnet::cv.glmnet(X_mat, y_vec, alpha = a, nfolds = 10)
  data.frame(alpha = a,
             lambda_min = cv$lambda.min,
             cvm_min    = min(cv$cvm))
})
en_grid <- do.call(rbind, en_results)
cat("\nElastic net: CV error by alpha\n")
round(en_grid, 4)

best_alpha <- en_grid$alpha[which.min(en_grid$cvm_min)]
best_lambda_en <- en_grid$lambda_min[which.min(en_grid$cvm_min)]
cat("Best alpha:", best_alpha, "| Best lambda:", round(best_lambda_en, 4), "\n")

en_fit<-glmnet::glmnet(X_mat, y_vec, alpha = best_alpha,lambda = best_lambda_en)
en_coef<-as.matrix(coef(en_fit))
cat("\nElastic net coefficients\n")
print(round(en_coef, 4))

ggplot(en_grid, aes(x = alpha, y = cvm_min)) +
  geom_line(linewidth = 1.1) + geom_point(size = 3) +
  labs(title = "Elastic net: CV MSE by alpha",
       x = "Alpha (0=ridge, 1=LASSO)", y = "CV MSE") +
  theme_bw(base_size = 14)

#- - - - - -
# Partial least square (PLS)

# Cross-validation to select number of components
pls_cv <- pls::plsr(y ~ ., data = sim_data, scale = TRUE,validation = "CV", segments = 10)
plot(pls_cv, "validation", estimate = "CV",
     main = "PLS: cross-validated RMSEP by number of components")

# Elbow / minimum CV error
cv_rmsep <- pls::RMSEP(pls_cv)$val[1, 1, ]
best_ncomp_pls <- which.min(cv_rmsep) - 1   # subtract 1 for intercept
cat("Optimal PLS components (CV):", best_ncomp_pls, "\n")

pls_fit <- pls::plsr(y ~ ., data = sim_data, scale = TRUE,
                     ncomp = best_ncomp_pls)
summary(pls_fit)

# Variance explained per component
cat("\nVariance explained (X and Y)\n")
pls::explvar(pls_fit)

# VIP scores (variable importance in projection)
vip_pls<-function(object) {
  SS<-c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm <-object$loading.weights /
    matrix(sqrt(colSums(object$loading.weights^2)),
           nrow = nrow(object$loading.weights), ncol = ncol(object$loading.weights),
           byrow = TRUE)
  SSW<- matrix(SS, nrow = nrow(Wnorm), ncol = length(SS), byrow = TRUE) * Wnorm^2
  sqrt(nrow(Wnorm) * rowSums(SSW) / sum(SS))
}
vip_scores<-vip_pls(pls_fit)
names(vip_scores)<-colnames(X_raw)

cat("\nPLS VIP scores (>1 = important)\n")
round(sort(vip_scores, decreasing = TRUE), 3)

ggplot(data.frame(var = names(vip_scores), vip = vip_scores),
       aes(x = reorder(var, vip), y = vip)) +
  geom_col(fill = "#A6DDCE", color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "PLS: Variable Importance in Projection (VIP)",
       x = "", y = "VIP score") +
  theme_bw(base_size = 13)


## ============================================================
## Comparison: cross-validated RMSE for all methods
## ============================================================
cat("\n\n===== MODEL COMPARISON: 10-fold CV RMSE =====\n")

k_folds  <- 10
fold_idx <- caret::createFolds(y_vec, k = k_folds, list = TRUE)

cv_rmse <- function(model_fn) {
  rmse_folds <- sapply(fold_idx, function(idx) {
    train <- sim_data[-idx, ]
    test  <- sim_data[idx, ]
    preds <- tryCatch(model_fn(train, test), error = function(e) NA)
    if (any(is.na(preds))) return(NA)
    sqrt(mean((test$y - preds)^2))
  })
  round(mean(rmse_folds, na.rm = TRUE), 4)
}

rmse_ols   <- cv_rmse(function(tr, te) predict(lm(y~., tr), te))
rmse_ridge <- cv_rmse(function(tr, te) {
  fit <- glmnet::glmnet(as.matrix(tr[,-1]), tr$y, alpha=0, lambda=lambda_ridge_1se)
  predict(fit, as.matrix(te[,-1]))[,1]
})
rmse_lasso <- cv_rmse(function(tr, te) {
  fit <- glmnet::glmnet(as.matrix(tr[,-1]), tr$y, alpha=1, lambda=lambda_lasso_1se)
  predict(fit, as.matrix(te[,-1]))[,1]
})
rmse_en <- cv_rmse(function(tr, te) {
  fit <- glmnet::glmnet(as.matrix(tr[,-1]), tr$y, alpha=best_alpha, lambda=best_lambda_en)
  predict(fit, as.matrix(te[,-1]))[,1]
})

comp_df <- data.frame(
  Method = c("OLS","Ridge","LASSO","Elastic Net"),
  CV_RMSE = c(rmse_ols, rmse_ridge, rmse_lasso, rmse_en)
)
cat("\n--- CV RMSE comparison ---\n")
print(comp_df)

ggplot(comp_df, aes(x = reorder(Method, CV_RMSE), y = CV_RMSE)) +
  geom_col(fill = "#A6DDCE", color = "black") +
  labs(title = "Model comparison: 10-fold CV RMSE",
       x = "", y = "CV RMSE") +
  theme_bw(base_size = 14)

#--------------------------------------------
## Reusable pipeline
run_linear_pipeline<-function(data, outcome_col, predictor_cols,
                                methods = c("ols","robust","lasso","ridge","enet","pls"),
                                nfolds  = 10, seed = 123) {
  set.seed(seed)
  df<-data[, c(outcome_col, predictor_cols)]
  X<-as.matrix(df[, predictor_cols])
  y<-df[[outcome_col]]
  form<-as.formula(paste(outcome_col, "~ ."))
  res<-list()
  
  if ("ols" %in% methods) res$ols<-lm(form, data = df)
  if ("robust"%in% methods) res$robust<- MASS::rlm(form, data = df, method="MM")
  if ("ridge" %in% methods) {
    cv_r<-glmnet::cv.glmnet(X, y, alpha=0, nfolds=nfolds)
    res$ridge <- glmnet::glmnet(X, y, alpha=0, lambda=cv_r$lambda.1se)
    res$ridge_lambda <- cv_r$lambda.1se
  }
  if ("lasso"  %in% methods) {
    cv_l<-glmnet::cv.glmnet(X, y, alpha=1, nfolds=nfolds)
    res$lasso<-glmnet::glmnet(X, y, alpha=1, lambda=cv_l$lambda.1se)
    res$lasso_lambda<-cv_l$lambda.1se
  }
  if ("enet"   %in% methods) {
    best_a<-sapply(seq(0.1,0.9,0.2), function(a) {
      min(glmnet::cv.glmnet(X, y, alpha=a, nfolds=nfolds)$cvm)
    })
    a_opt<-seq(0.1,0.9,0.2)[which.min(best_a)]
    cv_e<-glmnet::cv.glmnet(X, y, alpha=a_opt, nfolds=nfolds)
    res$enet<-glmnet::glmnet(X, y, alpha=a_opt, lambda=cv_e$lambda.1se)
  }
  if ("pls" %in% methods) {
    pls_cv_r<- pls::plsr(form, data=df, scale=TRUE, validation="CV", segments=nfolds)
    nc<-which.min(pls::RMSEP(pls_cv_r)$val[1,1,]) - 1
    nc<-max(nc, 1)
    res$pls<-pls::plsr(form, data=df, scale=TRUE, ncomp=nc)
  }
  
  cat("\nLinear pipeline complete. Methods fitted:", paste(names(res), collapse=", "), "\n")
  return(res)
}

all_models<-run_linear_pipeline(
  data=sim_data,
  outcome_col="y",
  predictor_cols= paste0("X", 1:p),
  methods=c("ols","robust","lasso","ridge","enet","pls"))

all_models
