# ------------------------------------------------------------
# Purpose: Load a saved model and use it to make predictions
# ------------------------------------------------------------

#---------------------------------------------------------
## 1. Loading the saved model

loaded_model <- readRDS("example_model.rds")

#---------------------------------------------------------
## 2. Preparing new data for prediction

# Replace this with your own new / unseen data
new_data <- data.frame(
  x1 = c(0.1, 0.5),
  x2 = c(-1.2, 0.3)
)

#---------------------------------------------------------
## 3. Making predictions 

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
