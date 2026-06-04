#-------------------------------------------------------------------------------
## Reproducible & Generalisable DAG Script
# - Builds and visualises DAGs
# - Identifies adjustment sets, backdoor paths, instruments
# - Tests DAG implications against data (d-separation tests)
# - Simulates data consistent with the DAG and demonstrates confounding
# - Includes a flexible build_dag() helper
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup
rm(list = ls())
set.seed(123)

required_pkgs <- c("dagitty", "ggdag", "ggplot2", "dplyr","MASS", "lavaan","performance")
is_installed<-required_pkgs %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed == FALSE)){
  install.packages(required_pkgs[!is_installed],repos = "http://cran.us.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

#--------------------------------------------
## Step 2: Defining DAG

# Example: environmental epidemiology
# Exposure (E) -> Outcome (Y)
# Confounder (C) -> E and Y
# Mediator (M): E -> M -> Y
# Instrument (Z): Z -> E only
# Collider (K): E -> K <- Y (do not adjust for K)

dag1<-dagitty::dagitty('dag {
  Z [pos="0,1"]
  C [pos="1,0"]
  E [pos="1,1" exposure]
  M [pos="2,1"]
  K [pos="2,2"]
  Y [pos="3,1" outcome]

  Z -> E
  C -> E
  C -> Y
  E -> M
  E -> K
  M -> Y
  E -> Y
  Y -> K
}')

#--------------------------------------------
## Step 3: Visualizing DAG

ggdag::ggdag(dag1, layout = "nicely") +
  ggdag::theme_dag() +
  labs(title = "DAG: exposure (E) -> outcome (Y) with confounding, mediation, collider")

# Highlighting paths
ggdag::ggdag_paths(dag1, from = "E", to = "Y") +
  ggdag::theme_dag() +
  labs(title = "All paths from E to Y")

#--------------------------------------------
## Step 4: Identifying adjustment sets

cat("\nMinimal adjustment sets (to estimate total E -> Y effect)\n")
dagitty::adjustmentSets(dag1, exposure = "E", outcome = "Y", type = "minimal")

cat("\nAll sufficient adjustment sets\n")
dagitty::adjustmentSets(dag1, exposure = "E", outcome = "Y", type = "all")

cat("\nAdjustment sets for DIRECT effect (E -> Y, blocking M)\n")
dagitty::adjustmentSets(dag1, exposure = "E", outcome = "Y",effect = "direct", type = "minimal")

#--------------------------------------------
## Step 5: Identifying instrumental variables

cat("\nInstrumental variables for E -> Y\n")
dagitty::instrumentalVariables(dag1, exposure = "E", outcome = "Y")

#--------------------------------------------
## Step 6: Identifying colliders and backdoor paths ----
cat("\nBackdoor paths from E to Y\n")
dagitty::paths(dag1, from = "E", to = "Y", directed = FALSE)

# Warning about collider K
cat("\nIs K a collider on path E -> K <- Y?\n")
cat("Conditioning on a collider OPENS a path — do not adjust for K.\n")

# Visualizing adjustment sets
ggdag::ggdag_adjustment_set(dag1, exposure = "E", outcome = "Y",shadow = TRUE) +
  ggdag::theme_dag() +
  labs(title = "Sufficient adjustment sets (shaded = must adjust)")

#--------------------------------------------
## Step 7: Simulating data

n<-1000
C<-rnorm(n)
Z<-rnorm(n)
E<-0.5*C + 0.6*Z + rnorm(n)
M<-0.4*E + rnorm(n)
Y<-0.3*E + 0.5*M + 0.4*C + rnorm(n)
K<-0.3*E + 0.3*Y + rnorm(n)   # collider

dag_data<-data.frame(Z=Z, C=C, E=E, M=M, K=K, Y=Y)

#--------------------------------------------
## Step 8: Demonstrating confounding and proper adjustment

# Unadjusted (biased — ignores C)
fit_unadj<-lm(Y ~ E, data = dag_data)
cat("\nUnadjusted E -> Y (biased)\n")
coef(fit_unadj)["E"]

# Adjusted for C (correct — blocks backdoor C -> E, C -> Y)
fit_adj <- lm(Y ~ E + C, data = dag_data)
cat("\nAdjusted for C (unbiased total effect)\n")
coef(fit_adj)["E"]

# check_dag(fit_adj,Y ~ E,outcome="Y",exposure="E",adjusted="C")
# plot(check_dag(fit_adj,Y ~ E,outcome="Y",exposure="E",adjusted="C"))
# ggdag::ggdag_status(check_dag(fit_adj,Y ~ E,outcome="Y",exposure="E",adjusted="C"))

# Over-adjusted: conditioning on mediator M (blocks part of effect)
fit_over <- lm(Y ~ E + C + M, data = dag_data)
cat("\n Over-adjusted (+ M): direct effect only\n")
coef(fit_over)["E"]

# Collider bias: conditioning on K (WRONG — opens E-Y path via K)
fit_collider <- lm(Y ~ E + C + K, data = dag_data)
cat("\nCollider bias (+ K, incorrect): biased estimate\n")
coef(fit_collider)["E"]

#--------------------------------------------
## Step 9: Testing DAG implications (d-separation tests)

# Conditional independence claims implied by the DAG
cat("\nImplied conditional independencies\n")
impl_ci <- dagitty::impliedConditionalIndependencies(dag1)
impl_ci

# Testing each implied CI against the data using partial correlations
# (works for linear / Gaussian data; use lavaan for SEM-based tests)
cat("\nTesting implied independencies against data\n")
test_res <- dagitty::localTests(dag1, data = dag_data, type = "cis")
round(test_res, 3)
cat("p < 0.05 = evidence against that independence claim (potential DAG mis-specification)\n")

#--------------------------------------------
## Step 10: Mediation analysis via DAG-guided approach
# Total effect = direct + indirect (through M)
# Using difference-in-coefficients method

total_eff<-coef(fit_adj)["E"]
direct_eff<-coef(fit_over)["E"]
indirect<-total_eff - direct_eff
cat("\nMediation decomposition\n")
cat("Total effect:   ", round(total_eff,  4), "\n")
cat("Direct effect:  ", round(direct_eff, 4), "\n")
cat("Indirect (via M):", round(indirect,  4), "\n")
cat("Proportion mediated:", round(indirect / total_eff * 100, 1), "%\n")

#--------------------------------------------
## Step 11: More complex DAG — time-varying confounding
dag_tv<-dagitty::dagitty('dag {
  L0 [pos="0,0"]
  A0 [pos="1,0" exposure]
  L1 [pos="2,0"]
  A1 [pos="3,0"]
  Y  [pos="4,0" outcome]

  L0 -> A0
  L0 -> L1
  A0 -> L1
  L1 -> A1
  A0 -> Y
  A1 -> Y
  L1 -> Y
}')

ggdag::ggdag(dag_tv, layout = "nicely") +
  ggdag::theme_dag() +
  labs(title = "Time-varying confounding DAG (g-methods needed)")

cat("\nAdjustment sets for A0+A1 -> Y (time-varying)\n")
dagitty::adjustmentSets(dag_tv, exposure = c("A0","A1"),outcome = "Y", type = "minimal")

#--------------------------------------------
## Step 12: Reusable DAG helper

build_dag<-function(dag_string, exposure, outcome,data = NULL, plot = TRUE){
  dag<-dagitty::dagitty(dag_string)
  dagitty::exposures(dag)<-exposure
  dagitty::outcomes(dag)<-outcome
  
  if(plot){
    p<-ggdag::ggdag(dag, layout = "nicely") +
      ggdag::theme_dag() +
      labs(title = paste("DAG:", exposure, "->", outcome))
    print(p)
    
    p2<-ggdag::ggdag_adjustment_set(dag, exposure = exposure,
                                      outcome = outcome, shadow = TRUE) +
      ggdag::theme_dag() +
      labs(title = "Adjustment sets")
    print(p2)
  }
  
  adj<-dagitty::adjustmentSets(dag, type = "minimal")
  inst<-tryCatch(dagitty::instrumentalVariables(dag, exposure, outcome),
                   error = function(e) NULL)
  
  result<-list(dag = dag, adjustment_sets = adj, instruments = inst)
  
  if(!is.null(data)){
    ci_tests<-dagitty::impliedConditionalIndependencies(dag)
    local_res<-tryCatch(dagitty::localTests(dag, data = data, type = "cis"),
                          error = function(e) NULL)
    result$local_tests<-local_res
  }
  
  return(result)
}

# Example
dag_res<-build_dag(dag_string = 'dag {
    C [pos="1,0"]; Z [pos="0,1"]
    E [pos="1,1" exposure]; M [pos="2,1"]
    Y [pos="3,1" outcome]; K [pos="2,2"]
    Z->E; C->E; C->Y; E->M; M->Y; E->Y; E->K; Y->K
  }',
  exposure="E",
  outcome="Y",
  data=dag_data)
