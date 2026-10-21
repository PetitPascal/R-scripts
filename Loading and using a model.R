# ------------------------------------------------------------
# Purpose: Load a saved model and use it to make predictions
# ------------------------------------------------------------

#----------------------------------------------------------------
#### Configurations ####

## Disabling memory torture
gctorture(FALSE)

## Installing and loading packages
pack_needed<-c("tidyverse","xgboost","qs2")

is_installed<-pack_needed %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed == FALSE)){
  install.packages(pack_needed[!is_installed],repos = "http://cran.us.r-project.org")
}
invisible(lapply(pack_needed, library, character.only = TRUE))

#----------------------------------------------------------------
#### 1. Loading the saved model ####

# RDS
loaded_model <- readRDS("Example_model.rds") # dependent of the package version, better solution below

# # QS2
# loaded_model <- qd_read("example_model.qs")

# # xgboost
# loaded_model <- xgb.load("example_model.ubj")

#----------------------------------------------------------------
#### 2. Preparing new data for prediction ####

# Replace this with your own new / unseen data
new_data <- data.frame(x1 = c(0.1, 0.5), x2 = c(-1.2, 0.3))

#----------------------------------------------------------------
#### 3. Making predictions ####

# Some models (e.g., xgboost) require a matrix, others (glm, randomForest, caret) require a data frame.
if (inherits(loaded_model, "xgb.Booster")) {
  # XGBoost requires a numeric matrix
  pred_input <- as.matrix(new_data)
} else {
  # Most other R models work with a data frame
  pred_input <- new_data
}

predictions <- predict(loaded_model, newdata = pred_input)
print(predictions)
