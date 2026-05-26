#-------------------------------------------------------------------------------
## Reproducible & Generalisable Hidden Markov Model (HMM) Script
# - Simulates longitudinal data with latent state switching
# - Fits HMM (continuous/categorical observations)
# - Tests assumptions: Viterbi decoding, state classification, BIC selection
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup ----
rm(list = ls())
set.seed(123)

required_pkgs<-c("depmixS4", "MASS", "dplyr", "ggplot2","tidyr", "RColorBrewer", "mclust")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

#--------------------------------------------
## Step 2: Simulating HMM data
# 2 latent states: "low" and "high" activity
# State-dependent observation: continuous variable (e.g. step count)
n_subj<-80
n_obs<-12   # observations per subject (e.g. monthly)
N<-n_subj * n_obs

# True transition matrix
A_true <-matrix(c(0.85, 0.15,   # from state 1
                  0.20, 0.80),  # from state 2
                  nrow = 2, byrow = TRUE)

# State-dependent emission means
mu_true<-c(3000, 8000)   # step counts
sd_true<-c(800,  1200)

# Simulating latent states and observations
states_all<-numeric(N)
obs_all<-numeric(N)
subj_all<-rep(seq_len(n_subj), each = n_obs)

for (s in seq_len(n_subj)) {
  idx<-((s-1)*n_obs + 1):(s*n_obs)
  st<-integer(n_obs)
  st[1]<-sample(1:2, 1, prob = c(0.6, 0.4))   # initial state
  for (t in 2:n_obs)
    st[t]<-sample(1:2, 1, prob = A_true[st[t-1], ])
  states_all[idx]<-st
  obs_all[idx]<-rnorm(n_obs, mu_true[st], sd_true[st])
}

sim_hmm<-data.frame(subject=factor(subj_all),
  time=rep(seq_len(n_obs), n_subj),
  y=pmax(obs_all, 0),    # non-negative
  true_state=states_all)

#--------------------------------------------
## Step 3: Assumption / EDA checks

# 3a. Distribution of observations
ggplot(sim_hmm, aes(x = y)) +
  geom_histogram(fill = "#A6DDCE", color = "black", bins = 40) +
  labs(title = "Observation distribution", x = "y", y = "Count") +
  theme_bw()

# 3b. Preliminary evidence for latent states using Gaussian mixture
mclust_fit<-mclust::Mclust(sim_hmm$y, G = 1:5)
cat("\nMclust BIC-selected number of components:", mclust_fit$G, "\n")
plot(mclust_fit, what = "BIC")

# 3c. Time series plot (first 5 subjects)
sim_hmm %>% filter(subject %in% levels(subject)[1:5]) %>%
  ggplot(aes(x = time, y = y, group = subject, color = subject)) +
  geom_line() + geom_point() +
  labs(title = "Observed trajectory (first 5 subjects)",
       x = "Time", y = "y") +
  theme_bw() + theme(legend.position = "bottom")

#--------------------------------------------
## Step 4: Fitting HMM
# ntimes: vector of observation lengths per subject
ntimes_vec<-rep(n_obs, n_subj)

# Model with 2 hidden states, Gaussian response
hmm_2<-depmixS4::depmix(response=y ~ 1,# response model: mean per state
  data=sim_hmm,
  nstates=2,
  ntimes=ntimes_vec,
  family=gaussian())

fit_2<-depmixS4::fit(hmm_2, verbose = FALSE)
summary(fit_2)

#--------------------------------------------
## Step 5: Model selection — BIC across different numbers of states
bic_vals<-sapply(1:4, function(k) {
  m<-depmixS4::depmix(y ~ 1, data = sim_hmm, nstates = k,
                        ntimes = ntimes_vec, family = gaussian())
  f<-tryCatch(depmixS4::fit(m, verbose = FALSE),
                error = function(e) NULL)
  if (is.null(f)) return(NA)
  BIC(f)
})

bic_df<-data.frame(n_states = 1:4, BIC = bic_vals)
bic_df

ggplot(bic_df, aes(x = n_states, y = BIC)) +
  geom_line(linewidth = 1.1) + geom_point(size = 3) +
  scale_x_continuous(breaks = 1:4) +
  labs(title = "HMM model selection: BIC by number of states",
       x = "Number of states", y = "BIC") +
  theme_bw(base_size = 14)

best_k<-bic_df$n_states[which.min(bic_df$BIC)]
cat("\nBIC-selected number of states:", best_k, "\n")

#--------------------------------------------
## Step 6: Fitting best model
hmm_best<-depmixS4::depmix(y ~ 1, data = sim_hmm, nstates = best_k,
                             ntimes = ntimes_vec, family = gaussian())
fit_best<-depmixS4::fit(hmm_best, verbose = FALSE)
summary(fit_best)

#--------------------------------------------
## Step 7: Viterbi decoding (i.e., most likely state sequence)
viterbi_states <-depmixS4::viterbi(fit_best)
sim_hmm$decoded<-viterbi_states$state

# Confusion with true states (if known)
if (best_k == 2) {
  cat("\n--- State agreement (Viterbi vs true states) ---\n")
  conf_mat<-table(True = sim_hmm$true_state, Decoded = sim_hmm$decoded)
  print(conf_mat)
  # Note: label switching may occur — states may be relabelled
  cat("Accuracy (best alignment):",
      round(max(sum(diag(conf_mat)), sum(anti_diag<-conf_mat[1,2] + conf_mat[2,1])) /
              sum(conf_mat), 3), "\n")
}

#--------------------------------------------
## Step 8: Posterior state probabilities ----
post_probs   <-depmixS4::posterior(fit_best)
sim_hmm$p_s1 <-post_probs[, 2]   # P(state 1)
sim_hmm$p_s2 <-if (best_k >= 2) post_probs[, 3] else NA

# Plotting decoded states for first 5 subjects
sim_hmm %>%
  filter(subject %in% levels(subject)[1:5]) %>%
  mutate(decoded = factor(decoded)) %>%
  ggplot(aes(x = time, y = y, color = decoded, group = subject)) +
  geom_line(aes(group = subject), color = "grey70") +
  geom_point(size = 2.5) +
  facet_wrap(~subject, nrow = 1) +
  labs(title = "Viterbi-decoded states (first 5 subjects)",
       x = "Time", y = "y", color = "State") +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 13)

#--------------------------------------------
## Step 9: Extracting transition and emission parameters
cat("\nTransition matrix\n")
# Extracting from fitted model
trans_mat<-matrix(getpars(fit_best)[seq_len(best_k^2)], nrow = best_k, byrow = TRUE)
rownames(trans_mat)<-colnames(trans_mat)<-paste0("state_", seq_len(best_k))
print(round(trans_mat, 3))

cat("\nEmission parameters (mean, SD per state)\n")
emit_pars<-getpars(fit_best)
# Parameters after transition block: intercepts and SDs for Gaussian responses
n_trans  <-best_k + best_k * best_k   # initial + transition params
emis_part<-emit_pars[(n_trans + 1):length(emit_pars)]
cat(round(emis_part, 3), "\n")

#--------------------------------------------
## Step 10: Assumption checks

# 10a. Residuals within state (should be roughly normal)
for (k in seq_len(best_k)) {
  r<-sim_hmm$y[sim_hmm$decoded == k] -
    mean(sim_hmm$y[sim_hmm$decoded == k])
  sw<-shapiro.test(r[seq_len(min(length(r), 5000))])
  cat("State", k, "residuals — Shapiro-Wilk p =", round(sw$p.value, 4), "\n")
}

# 10b. State prevalence — warn if any state < 5%
state_prev<-prop.table(table(sim_hmm$decoded))
cat("\nState prevalence\n"); print(round(state_prev, 3))
if (any(state_prev < 0.05))
  warning("At least one state has < 5% prevalence — consider fewer states.")

# 10c. Autocorrelation of residuals (should be low within decoded state)
acf(r, main = "ACF of HMM residuals", lag.max = 20)

#--------------------------------------------
## Step 11: Covariate-dependent transition probabilities
# Example: transition probabilities depend on age (time-constant covariate)
sim_hmm_cov<-sim_hmm
sim_hmm_cov$age_sc<-scale(rep(rnorm(n_subj, 50, 10), each = n_obs))

hmm_cov<-depmixS4::depmix(response= y ~ 1,
  transition=~ age_sc,   # covariate on transition probabilities
  data=sim_hmm_cov,
  nstates=best_k,
  ntimes=ntimes_vec,
  family=gaussian())
fit_cov<-tryCatch(depmixS4::fit(hmm_cov, verbose = FALSE),
                    error = function(e) {
                      message("Covariate HMM failed: ", e$message)
                      NULL
                    })
if (!is.null(fit_cov)) {
  cat("\n--- Covariate-dependent HMM ---\n")
  print(summary(fit_cov))
  cat("BIC (no covariate):", BIC(fit_best),
      " | BIC (with age):", BIC(fit_cov), "\n")
}

# #--------------------------------------------
# ## Step 12: Reusable pipeline
# run_hmm_pipeline<-function(data, outcome_col, subject_col, time_col,
#                              n_states_max = 4, family = gaussian(),
#                              covariate_cols = NULL, seed = 2025) {
#   set.seed(seed)
#   
#   # Sorting by subject then time
#   data<-data %>% arrange(.data[[subject_col]], .data[[time_col]])
#   
#   ntimes_v<-as.integer(table(data[[subject_col]]))
#   
#   # Model selection
#   bic_v<-sapply(seq_len(n_states_max), function(k) {
#     resp<-as.formula(paste(outcome_col, "~ 1"))
#     m   <-depmixS4::depmix(resp, data = data, nstates = k,
#                              ntimes = ntimes_v, family = family)
#     f<-tryCatch(depmixS4::fit(m, verbose = FALSE), error = function(e) NULL)
#     if (is.null(f)) return(NA)
#     BIC(f)
#   })
#   
#   best_k<-which.min(bic_v)
#   message("BIC-selected states: ", best_k)
#   
#   # Fitting best
#   resp_f<-as.formula(paste(outcome_col, "~ 1"))
#   m_best<-depmixS4::depmix(resp_f, data = data, nstates = best_k,
#                              ntimes = ntimes_v, family = family)
#   fit_b <-depmixS4::fit(m_best, verbose = FALSE)
#   
#   # Decode
#   data$decoded_state<-depmixS4::viterbi(fit_b)$state
#   data$posterior_s1 <-depmixS4::posterior(fit_b)[, 2]
#   
#   return(list(fit = fit_b, data = data, bic = bic_v, best_k = best_k))
# }
# 
# hmm_res<-run_hmm_pipeline(
#   data          = sim_hmm,
#   outcome_col   = "y",
#   subject_col   = "subject",
#   time_col      = "time",
#   n_states_max  = 4
# )
# cat("\nPipeline complete. Best k =", hmm_res$best_k, "\n")
