#----------------------------------------------------------------
#### Configurations ####

## Disabling memory torture
gctorture(FALSE)

## Installing packages if needed
pack_needed<-c("data.table","tidyverse","mllrnrs","broom","doParallel","foreach",
               "splitTools","conflicted","grid","gridExtra","RColorBrewer","mlbench",
               "mlexperiments","caret","MLmetrics","patchwork",
               "xgboost","parallel","here")
for (i in 1:length(pack_needed)){
  if(pack_needed[i]%in%.packages(all.available=TRUE)){
  }else{
    install.packages(pack_needed[i])
  }
}

## Package loading
library(tidyverse)
library(data.table)
library(broom)
library(doParallel)
library(foreach)
library(pROC)
library(grid)
library(gridExtra)
library(conflicted)
library(RColorBrewer)
library(mlexperiments)
library(mllrnrs)
library(mlbench)
library(caret)
library(MLmetrics)
library(xgboost)
library(patchwork)
library(here)

## Preventing package conflicts
conflict_prefer("select", "dplyr") 
conflict_prefer("filter", "dplyr") 
conflict_prefer("slice", "dplyr")
conflict_prefer("alpha", "scales")
conflict_prefer("auc", "pROC")
conflict_prefer("roc", "pROC")

## Setting the working directory
here::here("XGBoost - Binary classification - nested CV")

## Setting seed
seed <- 123

## Setting cores
if (isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_")))) {
  ncores <- 2L
} else {
  ncores <- ifelse(
    test = parallel::detectCores() > 4,
    yes = 4L,
    no = ifelse(
      test = parallel::detectCores() < 2L,
      yes = 1L,
      no = parallel::detectCores()
    )
  )
}

## Setting mlexperiments package options
options("mlexperiments.bayesian.max_init" = 10L)
options("mlexperiments.optim.xgb.nrounds" = 100L)
options("mlexperiments.optim.xgb.early_stopping_rounds" = 10L)

#----------------------------------------------------------------
#### Creating functions ####

## Function for plotting SHAP values
plot.shap.summary <- function(data_long){
  x_bound <- max(abs(data_long$value))
  require('ggforce')
  plot1 <- ggplot(data = data_long)+
    coord_flip() + 
    geom_sina(aes(x = variable, y = value, color = stdfvalue)) +
    geom_text(data = unique(data_long[, c("variable", "mean_value"), with = F]),
              aes(x = variable, y=-Inf, label = sprintf("%.3f", mean_value)),
              size = 3, alpha = 0.7,
              hjust = -0.2, 
              fontface = "bold") + # bold
    scale_color_gradient(low="#FFCC33", high="#6600CC", 
                         breaks=c(0,1), labels=c("Low","High")) +
    theme_bw() + 
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position="bottom") + 
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits = c(-x_bound, x_bound)) +
    scale_x_discrete(limits = rev(levels(data_long$variable)) 
    ) + 
    labs(y = "SHAP value (impact on model output)", x = "", color = "Feature value") 
  return(plot1)
}

## Function for standardizing feature values into the same range
std1 <- function(x){
  return ((x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T)))
}

## Function for formatting the display of summary statistics
test_format<-function(x){
  x<-as.numeric(x)
  sign_x<-if_else(x<0,"neg","pos")
  x<-abs(x)
  x_raw<-x
  
  if(is.na(x)|x==""|x=="NA"|is.infinite(x)){
    x_raw<-0
  }
  
  if(x_raw>=100){
    virg_pos<-str_locate(as.character(x_raw),"[.]")[1]
    
    if(!is.na(virg_pos)&as.numeric(substr(x_raw,virg_pos+1,virg_pos+1))>=5){
      x<-x+1
      x<-as.numeric(substr(x,1,virg_pos-1))
    }
    
  }
  
  if(is.na(x)|x==""|x=="NA"|is.infinite(x)){
    x<-""
  }else{
    
    if(x<0.01|x>=10000){
      x<-format(signif(x,3),scientific = T)
      
      if(nchar(x)==8&substr(x,4,4)!="e"){
        if(as.numeric(substr(x,4,4))>=5&x<0.01){
          
          x1<-as.numeric(substr(x,1,4))+0.1
          x2<-substr(x,5,8)
          x<-paste(substr(x1,1,4),x2,sep="")
        }else{
          x<-paste(substr(x,1,4),substr(x,5,8),sep="")
        }
      }else{
        x<-x
      }
      
    }else{
      
      x_save<-x
      x<-signif(x,3)
      
      if(nchar(x)==6){
        x<-as.numeric(substr(x,1,5))
      }
      
      if(nchar(x)==5){
        if(as.numeric(substr(x,5,5))>=5){
          x<-x+0.01
          x<-substr(x,1,4)
        }else{
          x<-substr(x,1,4)
        }
      }else{
        if(x>=1000){
          x<-as.character(signif(x_save,4))
        }else{
          x<-as.character(x)
        }
      }
    }
    
  }
  
  if(sign_x=="neg"&x_raw!=0){
    x<-paste("-",x,sep="")
  }
  
  if(x=="0e+00"){
    x<-"0"
  }
  
  return(x)
}

# Function for determining the direction (i.e.: how a feature impacts the model)
determine_direction <- function(df, feature_col, shap_col,
                                majority_threshold = 0.55,
                                n_bins = 200) {
  
  x <- df[[feature_col]]
  s <- df[[shap_col]]
  
  valid <- complete.cases(x, s)
  x <- x[valid]; s <- s[valid]
  
  if (length(unique(x)) <= 1 || length(unique(s)) <= 1)
    return("neutral")
  
  #---------------------------
  # BINARY / SPARSE FEATURES
  #---------------------------
  if (length(unique(x)) <= 2 || quantile(x, 0.75) == 0) {
    
    group0 <- s[x == min(x)]
    group1 <- s[x == max(x)]
    
    n0 <- length(group0); n1 <- length(group1)
    if (n0 == 0 || n1 == 0) return("neutral")
    
    group0_sorted <- sort(group0)
    
    pos_count <- sum(findInterval(group1, group0_sorted, left.open = TRUE))
    n_total <- as.double(n0) * as.double(n1)
    
    prop <- pos_count / n_total
    if (prop >= majority_threshold) return("promoting predictor")
    if (prop <= (1 - majority_threshold)) return("mitigating predictor")
    return("neutral")
  }
  
  #---------------------------
  # CONTINUOUS FEATURES
  #---------------------------
  
  bins <- cut(x, breaks = n_bins, include.lowest = TRUE)
  
  bx <- tapply(x, bins, mean)
  bs <- tapply(s, bins, mean)
  bc <- as.numeric(tapply(s, bins, length))
  
  keep <- !is.na(bx)
  bx <- bx[keep]
  bs <- bs[keep]
  bc <- bc[keep]
  
  B <- length(bx)
  if (B < 2) return("neutral")
  
  # Vectorized pairwise comparisons
  bx_mat_i <- matrix(bx, B, B)
  bx_mat_j <- t(bx_mat_i)
  
  bs_mat_i <- matrix(bs, B, B)
  bs_mat_j <- t(bs_mat_i)
  
  # Determining which bin has higher feature value
  higher_x <- bx_mat_i > bx_mat_j
  
  # Counting sample pairs
  bc_mat_i <- matrix(as.numeric(bc), B, B)
  bc_mat_j <- t(bc_mat_i)
  n_pairs  <- bc_mat_i * bc_mat_j
  
  total_pairs <- sum(n_pairs[higher_x])
  
  if (total_pairs == 0) return("neutral")
  
  pos_pairs <- sum(n_pairs[higher_x & (bs_mat_i > bs_mat_j)])
  neg_pairs <- sum(n_pairs[higher_x & (bs_mat_i < bs_mat_j)])
  
  prop_pos <- pos_pairs / total_pairs
  prop_neg <- neg_pairs / total_pairs
  
  if (prop_pos >= majority_threshold) return("promoting predictor")
  if (prop_neg >= majority_threshold) return("mitigating predictor")
  return("neutral")
}

## not including function (opposite function of %in%)
`%ni%`<-Negate('%in%')

#----------------------------------------------------------------
#### Data import and preprocessing ####

## Importing the clean dataset
clean_data<-as_tibble(na.omit(survival::colon) %>% filter(etype == 2) %>% select(-c(id,time,study,etype)) %>% mutate(rx=case_when(rx=="Obs"~0,rx=="Lev"~1,T~2)))

## Renaming the outcome variable
colnames(clean_data)[which(colnames(clean_data)=="status")]<-"outcome"
clean_data_save<-clean_data

## Creating a table with the variable type
var_type_tab<-as_tibble(data.frame(feature=colnames(clean_data),
                                   Var_type=c("ordered","dichotomous","continuous","dichotomous","dichotomous","dichotomous","continuous","dichotomous",
                                              "ordered","ordered","dichotomous","dichotomous")))

## Reorganising the dataset to ensure the outcome is the first column
clean_data<-clean_data %>% select(outcome,colnames(clean_data)[which(colnames(clean_data) %ni% c("outcome"))])

#--------------------------
## Transforming the clean dataset into data.table
dataset <- clean_data |> data.table::as.data.table()

## Creating vectors with the independent and dependent variable names
feature_cols <- colnames(dataset)[2:ncol(dataset)]
target_col <- "outcome"

## Identifying unordered factors
categorical_var<-var_type_tab %>% filter(feature %in% colnames(clean_data)) %>% 
  filter(Var_type %in% c("categorical")) %>% select(feature) %>% pull

## One-hot encoding of unordered factors
if (length(categorical_var) > 0) {
  for (col in categorical_var) {
    dataset[[col]] <- addNA(dataset[[col]])
  }
  dummies <- model.matrix(~ -1 + ., data = dataset[, ..categorical_var])
  dataset <- cbind(dataset[, !..categorical_var], dummies)
}

## Preparing training and test datasets
X_mat <- as.matrix(dataset[, feature_cols, with = FALSE])
y_vec <- dataset[[target_col]]


#----------------------------------------------------------------
#### 50x repeated 5-fold nested cross-validation ####

## Settings
k_folds_outer <- 5    # outer CV
k_folds_inner <- 3    # inner CV
n_repeats <- 50       # number of repeats
seed <- 1234          # random seed
all_results <- list() # for saving all results
counter <- 1

#--------------------------
## Performing the nested CV

set.seed(seed)
for (r in 1:n_repeats) {
  
  # Outer folds (for unbiased evaluation)
  outer_folds <- splitTools::create_folds(y = y_vec, k = k_folds_outer, type = "stratified", seed = seed + r)
  
  for (fold_idx in seq_along(outer_folds)) {
    
    test_idx <- outer_folds[[fold_idx]]                     # outer test set
    train_idx <- setdiff(seq_len(nrow(dataset)), test_idx)  # outer train set
    
    X_outer_train <- X_mat[train_idx, , drop = FALSE]
    y_outer_train <- y_vec[train_idx]
    X_outer_test  <- X_mat[test_idx, , drop = FALSE]
    y_outer_test  <- y_vec[test_idx]
    
    # Handling class imbalance in outer train
    weight <- max(table(y_outer_train)) / min(table(y_outer_train))
    
    # Inner CV for hyperparameter tuning
    inner_folds <- splitTools::create_folds(y = y_outer_train, k = k_folds_inner, type = "stratified", seed = seed + r + fold_idx)
    
    # Hyperparameter grid
    param_list <- expand.grid(
      subsample = seq(0.5, 1, 0.25),
      colsample_bytree = seq(0.5, 1, 0.25),
      min_child_weight = c(1, 5, 10),
      learning_rate = c(0.05, 0.01, 0.3),
      max_depth = c(3, 5, 7),
      scale_pos_weight = weight) %>% dplyr::slice_sample(n = 30)
    
    best_auc <- -Inf
    best_params <- NULL
    
    for (i in 1:nrow(param_list)) {
      aucs_inner <- c()
      
      for (inner_idx in seq_along(inner_folds)) {
        val_idx <- inner_folds[[inner_idx]]
        train_idx_inner <- setdiff(seq_len(nrow(X_outer_train)), val_idx)
        
        X_inner_tr <- X_outer_train[train_idx_inner, , drop = FALSE]
        y_inner_tr <- y_outer_train[train_idx_inner]
        X_inner_val <- X_outer_train[val_idx, , drop = FALSE]
        y_inner_val <- y_outer_train[val_idx]
        
        model_inner <- xgboost::xgb.train(
          params = as.list(param_list[i, ]),
          data = xgboost::xgb.DMatrix(X_inner_tr, label = y_inner_tr),
          nrounds = 100,
          objective = "binary:logistic",
          eval_metric = "auc",
          verbose = 0
        )
        
        pred_inner <- predict(model_inner, X_inner_val)
        auc_inner <- as.numeric(pROC::auc(y_inner_val, pred_inner))
        aucs_inner <- c(aucs_inner, auc_inner)
      }
      
      mean_inner_auc <- mean(aucs_inner, na.rm = TRUE)
      if (!is.na(mean_inner_auc) && mean_inner_auc > best_auc) {
        best_auc <- mean_inner_auc
        best_params <- param_list[i, ]
      }
    }
    
    # Training final model on full outer training set using best hyperparameters
    final_model <- xgboost::xgb.train(
      params = as.list(best_params),
      data = xgboost::xgb.DMatrix(X_outer_train, label = y_outer_train),
      nrounds = 100,
      objective = "binary:logistic",
      eval_metric = "auc",
      verbose = 0
    )
    
    # Predicting on outer test set
    pred_prob <- predict(final_model, X_outer_test)
    pred_bin  <- ifelse(pred_prob > 0.5, 1, 0)
    
    # Computing performance metrics
    fold_metrics <- data.frame(
      rep = r,
      fold = fold_idx,
      AUROC = as.numeric(pROC::auc(y_outer_test, pred_prob)),
      Accuracy = mean(pred_bin == y_outer_test),
      Sensitivity = MLmetrics::Recall(pred_bin, y_outer_test),
      Specificity = MLmetrics::Specificity(pred_bin, y_outer_test),
      F1 = MLmetrics::F1_Score(pred_bin, y_outer_test),
      PPV = MLmetrics::Precision(pred_bin, y_outer_test),
      NPV = sum((pred_bin == 0 & y_outer_test == 0))/sum(pred_bin == 0)
    )
    
    all_results[[counter]] <- fold_metrics
    counter <- counter + 1
  }
}

# Combining results into a summary table
all_results_df <- do.call(rbind, all_results)
as_tibble(all_results_df)

# Aggregating metrics to get a summary
metrics <- c("AUROC", "Accuracy", "Sensitivity", "Specificity", "F1", "PPV", "NPV")
summary_metrics <- lapply(metrics, function(m) {
  x <- all_results_df[[m]]
  tibble(
    Metric = m,
    Mean = mean(x, na.rm = TRUE),
    SD = sd(x, na.rm = TRUE),
    CI_lower = quantile(x, 0.025, na.rm = TRUE),
    CI_upper = quantile(x, 0.975, na.rm = TRUE))
}) %>% bind_rows()

# Brier Score: measures mean squared error of predicted probabilities
brier <- mean((pred_prob - y_outer_test)^2)
brier

#--------------------------
## Saving XGBoost results
write.table(all_results_df,paste("All CV results_",Sys.Date(),".csv"),sep=";",col.names=T,row.names=F)
write.table(summary_metrics,paste("Summary metrics_",Sys.Date(),".csv"),sep=";",col.names=T,row.names=F)
write.table(as_tibble(brier),paste("Brier_",Sys.Date(),".csv"),sep=";",col.names=T,row.names=F)

#----------------------------------------------------------------
#### SHAP analysis ####

#- - - - - - - - - -
## Looping over repeats and folds

shap_list <- list()
counter <- 1

for (r in 1:n_repeats) {
  folds <- splitTools::create_folds(y = y_vec, k = k_folds, type = "stratified", seed = seed + r)
  
  for (fold_idx in seq_along(folds)) {
    test_idx <- folds[[fold_idx]]
    train_idx <- setdiff(seq_len(nrow(dataset)), test_idx)
    
    # Splitting data
    X_tr <- X_mat[train_idx, , drop = FALSE]
    y_tr <- y_vec[train_idx]
    X_te <- X_mat[test_idx, , drop = FALSE]
    
    # Training model on this fold with tuned hyperparameters
    final_model <- xgboost::xgb.train(
      params = as.list(best_params),
      data = xgboost::xgb.DMatrix(X_tr, label = y_tr),
      nrounds = 100,
      objective = "binary:logistic",
      eval_metric = "auc",
      verbose = 0
    )
    
    # Computing SHAP values on validation set
    shap_values <- predict(final_model, X_te, predcontrib = TRUE)
    
    # Removing bias column
    shap_values <- shap_values[, !grepl("bias", colnames(shap_values), ignore.case = TRUE)]
    
    # Storing SHAP values in long format
    shap_dt <- as.data.table(shap_values)
    shap_dt[, sample_id := test_idx]
    shap_long <- melt(shap_dt, id.vars = "sample_id", variable.name = "feature", value.name = "shap_value")
    
    # Matching feature values
    fv_dt <- as.data.table(X_te)
    fv_dt[, sample_id := test_idx]
    fv_long <- melt(fv_dt, id.vars = "sample_id", variable.name = "feature", value.name = "feature_value")
    
    shap_long <- left_join(shap_long, fv_long, by = c("sample_id", "feature"))
    
    # Standardizing feature values for plotting
    shap_long <- shap_long %>% group_by(feature) %>%
      mutate(stdfvalue = (feature_value - min(feature_value, na.rm = TRUE)) /
               (max(feature_value, na.rm = TRUE) - min(feature_value, na.rm = TRUE))) %>%
      ungroup()
    
    shap_list[[counter]] <- shap_long
    counter <- counter + 1
  }
}

#- - - - - - - - - -
## Combining all folds and repeats
shap_all <- rbindlist(shap_list)

#- - - - - - - - - -
## Computing summary statistics per feature
shap_summary <- shap_all %>%
  group_by(feature) %>%
  summarise(
    mean_abs_shap = mean(abs(shap_value), na.rm = TRUE),
    sd_abs_shap = sd(abs(shap_value), na.rm = TRUE),
    CI_2.5 = quantile(abs(shap_value), 0.025, na.rm = TRUE),
    CI_97.5 = quantile(abs(shap_value), 0.975, na.rm = TRUE),
    mean_shap_signed = mean(shap_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_shap))

#- - - - - - - - - -
## Computing SHAP directionality
direction_impact <- shap_all %>%
  group_by(feature) %>%
  group_map(~ tibble(
    feature = .y$feature,
    direction = determine_direction(.x, "feature_value", "shap_value")
  )) %>%
  bind_rows()

shap_summary <- left_join(shap_summary, direction_impact, by = "feature")
shap_summary_save<-shap_summary
shap_summary<-shap_summary %>% filter(mean_abs_shap!=0) %>% 
  mutate(mean_abs_shap=if_else(direction=="limiting factor",-mean_abs_shap,mean_abs_shap),
         CI_2.5=if_else(direction=="limiting factor",-CI_2.5,CI_2.5),
         CI_97.5=if_else(direction=="limiting factor",-CI_97.5,CI_97.5))

#- - - - - - - - - -
## Creating long-format tibble for plotting
shap_plot_df <- shap_all %>%
  left_join(direction_impact, by = "feature") %>%
  mutate(impact = direction)

#- - - - - - - - - -
## Reordering features for plotting
ordre_def1 <- shap_summary %>% filter(mean_abs_shap >= 0) %>% arrange(desc(CI_2.5))
ordre_def2 <- shap_summary %>% filter(mean_abs_shap < 0) %>% arrange(desc(CI_97.5))
ordre_def <- bind_rows(ordre_def1, ordre_def2)
shap_summary$feature <- factor(shap_summary$feature, levels = rev(unique(ordre_def$feature)))
shap_summary<-shap_summary %>% mutate(pos=if_else(CI_97.5>=0,0.15,-0.15))
shap_summary<-shap_summary %>% rowwise %>% mutate(label_shap=test_format(abs(CI_97.5)),
                                                  label_mean_shap=test_format(abs(mean_abs_shap))) %>% ungroup

#- - - - - - - - - -
## Plot 1: error bar
plot_errorbar <- ggplot(shap_summary, aes(y = feature, x = mean_abs_shap, color = direction)) +
  geom_point() +
  geom_errorbar(aes(xmin = CI_2.5, xmax = CI_97.5), width = 0.2) +
  scale_x_continuous("mean |SHAP value|", labels = function(x) abs(x))+
  geom_text(aes(x=CI_97.5+pos,y=feature,label=label_shap),size=5)+
  scale_color_manual("Direction:",
                     values = c("promoting predictor" = "#C35C33",
                                "mitigating predictor" = "#40B696",
                                "neutral" = "black")) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0),
        strip.background = element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y = element_text(size=16,color="black"),
        axis.text.x = element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.position="bottom",
        axis.line = element_line(color = "black",size = 0.1, linetype = "solid"))

# Export
ggsave(plot_errorbar, file = paste0("SHAP_errorbar_", Sys.Date(), ".pdf"),
       dpi = 600, width = 60, height = 30, units = "cm", limitsize = FALSE)

#- - - - - - - - - -
## Plot 2: barplot 
shap_summary <- shap_summary %>% arrange(desc(mean_abs_shap))
shap_summary$feature <- factor(shap_summary$feature, levels = rev(shap_summary$feature))

plot_bar <- ggplot(shap_summary, aes(y = feature, x = mean_abs_shap, fill = direction)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(x=mean_abs_shap,y=feature,label=label_mean_shap),
            size=6,position=position_nudge(x = if_else(shap_summary$mean_abs_shap>=0,0.2,-0.2)))+
  scale_x_continuous("mean |SHAP value|", labels = function(x) abs(x)) +
  scale_fill_manual("Direction:",
                    values = c("mitigating predictor" = "#A6DDCE",
                               "promoting predictor" = "#F9CBC2",
                               "neutral" = "white")) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0),
        strip.background = element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y = element_text(size=16,color="black"),
        axis.text.x = element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.position="bottom",
        axis.line = element_line(color = "black",size = 0.1, linetype = "solid"))

# Export
ggsave(plot_bar, file = paste0("SHAP_barplot_", Sys.Date(), ".pdf"),
       dpi = 600, width = 60, height = 30, units = "cm", limitsize = FALSE)

#----------------------------------------------------------------
#### Saving SHAP results

# SHAP values
fwrite(shap_summary_save,paste("SHAP summary_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE)
fwrite(shap_all,paste("shap_all_",Sys.Date(),".txt",sep=""), sep = ";", row.names=FALSE)
