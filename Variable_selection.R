#-------------------------------------------------------------------------------
## Reproducible & Generalisable Variable Selection Script
# Covers: stepAIC, fastbw (rms), ols_step_best_subset (olsrr),GUESS, DSA algorithm, VarSelLCM 
# - Simulates data with sparse true signal
# - Runs each method and compares selected variable sets
# - Tests stability via bootstrap
#-------------------------------------------------------------------------------

#--------------------------------------------
## Setup
rm(list = ls())
set.seed(123)

required_pkgs <- c("MASS","dplyr","ggplot2","tidyr","caret","rms","olsrr","VarSelLCM","BeSS","glmnet","mombf")

is_installed<-required_pkgs %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed == FALSE)){
  install.packages(required_pkgs[!is_installed],repos = "http://cran.us.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

if (!requireNamespace("partDSA", quietly=TRUE))
  tryCatch(install.packages("partDSA"),
           error=function(e) message("partDSA not available on CRAN — skip DSA section"))

if(!("sparseMatrixStats"%in%.packages(all.available=TRUE))){
  BiocManager::install("sparseMatrixStats")
}

#--------------------------------------------
## Simulating data
n<-400
p_vars<-20

rho<-0.4
Sigma<-rho ^ as.matrix(dist(1:p_vars))
X<- as.data.frame(MASS::mvrnorm(n, mu=rep(0,p_vars), Sigma=Sigma))
colnames(X)<-paste0("V", seq_len(p_vars))

# True signal: only V1, V2, V3, V5, V8
beta_true<-rep(0, p_vars)
beta_true[c(1,2,3,5,8)]<-c(1.5,-1.0,0.8,0.6,-0.7)
y_cont<-as.numeric(as.matrix(X) %*% beta_true + rnorm(n, 0, 2))
y_bin<-rbinom(n, 1, plogis(as.matrix(X) %*% beta_true * 0.5))

sim_data<-cbind(data.frame(y=y_cont, y_bin=y_bin), X)
true_vars<-c("V1","V2","V3","V5","V8")

cat("True signal variables:", paste(true_vars, collapse=", "), "\n")
cat("Outcome R2 (approx):",round(cor(y_cont, as.matrix(X) %*% beta_true)^2, 3), "\n")

## Helper: evaluate selection
eval_selection <- function(selected, true_vars, all_vars) {
  tp<-sum(selected %in% true_vars)
  fp<-sum(!(selected %in% true_vars))
  fn<-sum(!(true_vars %in% selected))
  pre<-if(length(selected)>0) tp/(tp+fp) else 0
  rec<-tp / (tp + fn)
  f1<-if(pre+rec>0) 2*pre*rec/(pre+rec) else 0
  cat("Selected:", paste(sort(selected), collapse=", "), "\n")
  cat("TP:", tp, "| FP:", fp, "| FN:", fn,
      "| Precision:", round(pre,3), "| Recall:", round(rec,3),
      "| F1:", round(f1,3), "\n")
  invisible(data.frame(tp=tp, fp=fp, fn=fn,
                       precision=pre, recall=rec, f1=f1))
}

all_selected <- list()

#--------------------------------------------
## Method 1 - stepAIC (bidirectional)

full_model<-lm(y ~ ., data=sim_data[, c("y", paste0("V",1:p_vars))])
null_model<-lm(y ~ 1, data=sim_data)

step_aic<-MASS::stepAIC(full_model, direction="both",
                          scope=list(lower=null_model, upper=full_model),
                          trace=FALSE)

selected_aic<-names(coef(step_aic))[-1] # remove intercept
cat("stepAIC selected variables:\n")
all_selected$stepAIC <- eval_selection(selected_aic, true_vars, paste0("V",1:p_vars))

#--------------------------------------------
## Method 2 - fastbw (uses backward elimination on rms OLS/lrm models)

dd<-rms::datadist(sim_data)
options(datadist = "dd")

rms_ols<-rms::ols(y ~ ., data=sim_data[, c("y", paste0("V",1:p_vars))],x=TRUE, y=TRUE)

fbw_res<-rms::fastbw(rms_ols, rule="aic", type="residual")

selected_fbw<-fbw_res$names.kept
all_selected$fastbw<-eval_selection(selected_fbw, true_vars, paste0("V",1:p_vars))

#--------------------------------------------
## Method 3 — ols_step_best_subset (exhaustive search which is computationally expensive for p > ~15)

sub_data<-sim_data[, c("y", paste0("V", 1:12))]
ols_full<-lm(y ~ ., data=sub_data)

best_sub<-olsrr::ols_step_best_subset(ols_full)

# Selecting model with lowest BIC
best_k<-which.min(best_sub$metrics$sbc)
selected_bs<-strsplit(best_sub$metrics$predictors[best_k], " ")[[1]]
cat("Best subset (min BIC, k=", best_k, ") variables:\n")
all_selected$best_subset<-eval_selection(selected_bs, true_vars, paste0("V",1:12))

# Plot BIC vs model size
plot(best_sub$metrics$n, best_sub$metrics$sbc, type="b",
     xlab="Number of predictors", ylab="BIC",
     main="Best subset: BIC by model size")
abline(v=best_k, lty=2, col="red")

#--------------------------------------------
## Method 4 - LASSO (glmnet)

X_mat<-as.matrix(sim_data[, paste0("V",1:p_vars)])
cv_lasso<-glmnet::cv.glmnet(X_mat, y_cont, alpha=1, nfolds=10)
lasso_coef<-coef(cv_lasso, s="lambda.1se")
selected_lasso<-rownames(lasso_coef)[lasso_coef[,1] != 0 &rownames(lasso_coef) != "(Intercept)"]
all_selected$LASSO<-eval_selection(selected_lasso, true_vars, paste0("V",1:p_vars))

#--------------------------------------------
## Method 5 — VarSelLCM (for mixture/clustering-based selection): it identifies variables that discriminate latent classes

if (requireNamespace("VarSelLCM", quietly=TRUE)) {
  varsel_res <- tryCatch({
    VarSelLCM::VarSelCluster(
      x    = as.data.frame(X_mat[, 1:10]),   # subset for speed
      vbleSelec = TRUE,
      gvals=3,
      crit.varsel = "BIC")
  }, error=function(e){
    message("VarSelLCM failed: ", e$message); NULL
  })
  
  if(!is.null(varsel_res)) {
    relevant_vars <- varsel_res@model@names.relevant
    cat("VarSelLCM relevant variables:",
        paste(relevant_vars, collapse=", "), "\n")
  }
} else {
  message("VarSelLCM not available — skipping")
}

#--------------------------------------------
## Method 6 — DSA (deletion/substitution/addition)

  dsa_res<-tryCatch({
    partDSA::partDSA(x=X_mat[, 1:12],
      y=y_cont,
      control=partDSA::DSA.control(
        MPD=0.1,
        vfold=5)
    )
  }, error=function(e){
    message("partDSA failed: ", e$message); NULL
  })
  
  if (!is.null(dsa_res)) {
    cat("\nDSA results\n")
    print(dsa_res)
  }

#--------------------------------------------
## Method 7 - Graphical Unit Evolutionary Stochastic Search (GUESS) (Bayesian model selection with non-local priors) 

ms_res <- tryCatch({
    mombf::modelSelection(
      y     = y_cont,
      x     = X_mat[, 1:12],
      priorCoef  = mombf::momprior(tau=0.348),
      priorDelta = mombf::modelbbprior(alpha.p=1, beta.p=1),
      niter = 5000,
      verbose = FALSE
    )
  }, error=function(e){
    message("mombf failed: ", e$message); NULL
  })
  
  if (!is.null(ms_res)) {
    cat("\nmarginal inclusion probabilities\n")
    mip <- ms_res$margpp
    print(round(mip, 3))
    
    selected_mombf <- names(mip)[mip > 0.5]
    cat("Selected (MIP > 0.5):", paste(selected_mombf, collapse=", "), "\n")
    all_selected$mombf <- eval_selection(
      selected_mombf, true_vars, paste0("V",1:12)
    )
  }

#--------------------------------------------
## Bootstrap stability of variable selection - stepAIC

n_boot<-100
boot_select<-replicate(n_boot, {
  idx<-sample(seq_len(n), n, replace=TRUE)
  d_b<-sim_data[idx, c("y", paste0("V",1:p_vars))]
  full<-lm(y~., data=d_b)
  null<-lm(y~1, data=d_b)
  s<-tryCatch(
    names(coef(MASS::stepAIC(full, direction="both",
                             scope=list(lower=null,upper=full),
                             trace=FALSE)))[-1],
    error=function(e) character(0))
  s
}, simplify=FALSE)

select_freq<-sort(table(unlist(boot_select)) / n_boot,decreasing=TRUE)
cat("\nBootstrap inclusion frequencies (stepAIC)\n")
print(round(select_freq, 3))

freq_df<-data.frame(
  variable=names(select_freq),
  frequency=as.numeric(select_freq),
  true_var=names(select_freq) %in% true_vars)

ggplot(freq_df, aes(x=reorder(variable,frequency), y=frequency,
                    fill=true_var)) +
  geom_col(color="black") +
  geom_hline(yintercept=0.5, linetype="dashed", color="red") +
  coord_flip() +
  scale_fill_manual(values=c("TRUE"="#A6DDCE","FALSE"="#d6604d"),
                    labels=c("TRUE"="True signal","FALSE"="Noise")) +
  labs(title="Bootstrap inclusion frequency: stepAIC",
       subtitle="Red dashed = 50% threshold",
       x="", y="Inclusion frequency", fill="") +
  theme_bw(base_size=13)

#--------------------------------------------
## Summary comparison across methods
cat("True signal variables:", paste(true_vars, collapse=", "), "\n\n")

result_summary <- do.call(rbind,
                          lapply(names(all_selected), function(nm) {
                            cbind(method=nm, all_selected[[nm]])
                          })
)
print(result_summary)
