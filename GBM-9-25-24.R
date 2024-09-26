# GBM Model

library(gbm)
library(survival)
library(survcomp)
library(pdp)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(mice)

setwd("C:/Users/rmehta/OneDrive - Fred Hutchinson Cancer Research Center/BoxMigration/0. Fred Hutch Research/CIBMTR-Haplo-LEADER-PFS-Machine/1. RSM Data")

# Load the haven package
library(haven)


# Load the data from the CSV file
data <- read.csv("data_export.csv")

# Check data structure to confirm tpfs and pfs
str(data) 

# Check the first few rows of the data
head(data)



# Define your features and target
features <- c("ragecat_rsm", "dnrage", "rsex", "donorcmv", "recipientcmv", 
              "relative_rsm", "hctcigp", "b_leader", "drb1_gvh", "dqb1_gvh", "tce", "disease_stage_rsm")
target <- c("tpfs", "pfs") 

data$pfs <- data$pfs == 1 # Will be TRUE for events, FALSE otherwise

# Prepare the data for survival analysis
data$SurvObj <- Surv(data$tpfs, data$pfs == 1)

data$ragecat_rsm <- factor(data$ragecat_rsm)
data$rsex <- factor(data$rsex)
data$donorcmv <- factor(data$donorcmv)
data$recipientcmv <- factor(data$recipientcmv)
data$relative_rsm <- factor(data$relative_rsm)
data$hctcigp <- factor(data$hctcigp)
data$b_leader <- factor(data$b_leader)
data$drb1_gvh <- factor(data$drb1_gvh)
data$dqb1_gvh <- factor(data$dqb1_gvh)
data$tce <- factor(data$tce)
data$disease_stage_rsm <- factor(data$disease_stage_rsm)  # Include "disease_stage_rsm" as a factor

# REMOVE CASES WITH MISSING PFS/ TPFS INFO
data_gbm <- data[complete.cases(data[, c("tpfs", "pfs")]), ] 

# Fit the GBM model with cross-validation
gbm_model <- gbm(SurvObj ~ ., 
                 data = data_gbm[, c(features, "SurvObj")], 
                 distribution = "coxph", 
                 n.trees = 1000, 
                 interaction.depth = 3, 
                 cv.folds = 5)

# Get the optimal number of trees
n_trees <- gbm.perf(gbm_model, method = "cv")

print(n_trees)

# Print the model summary
summary(gbm_model)

# Variable importance
summary(gbm_model, cBars = 15, method = relative.influence, las = 2)

# C-index (Use n_trees in prediction)
risk_scores <- predict(gbm_model, newdata = data_gbm[, features], n.trees = n_trees)

c_index <- concordance.index(x = risk_scores, 
                             surv.time = data_gbm$tpfs, 
                             surv.event = data_gbm$pfs == 1)
print(c_index$c.index)

######################################################################
# Partial Dependence Plots
pdp_dnrage <- partial(gbm_model, pred.var = "dnrage", n.trees = n_trees, plot = TRUE)
plot(pdp_dnrage)

# Example: ICE Plot for 'dnrage'
ice_dnrage <- gbm_model %>%
  partial(pred.var = "dnrage", n.trees = n_trees, grid.resolution = 100, ice = TRUE) %>%
  autoplot(rug = TRUE, train = data, alpha = 0.1) +
  ggtitle("ICE Plot for Donor Age") +
  scale_y_continuous(labels = scales::comma)
ice_dnrage

######################################################################


# Example: LIME for a few observations
local_obs <- data[1:2, features]
explainer <- lime::lime(data[, features], gbm_model, 
                        model_type = "regression", 
                        predict_model = function(x, newdata){predict(x, newdata, n.trees = n_trees)})
explanation <- lime::explain(local_obs, explainer, n_features = 10) 
plot_explanations(explanation)


######################################################################
######################################################################
######################################################################

# Define the grid of hyperparameters to search
gbm_grid <- expand.grid(
  interaction.depth = c(2, 3, 4, 5),       
  n.trees = c(50, 100, 200, 500, 1000), 
  shrinkage = c(0.01, 0.05, 0.1),           
  n.minobsinnode = c(10, 20, 30)           
)

# Set up cross-validation
set.seed(123) 
cv_folds <- 5    
folds <- sample(1:cv_folds, size = nrow(data_gbm), replace = TRUE)

# Store results
results <- data.frame(interaction.depth = integer(),
                      n.trees = integer(),
                      shrinkage = double(),
                      n.minobsinnode = integer(),
                      c_index = double())

# Grid search loop 
for (i in 1:nrow(gbm_grid)) {
  # Get hyperparameters for this iteration
  depth <- gbm_grid$interaction.depth[i]
  trees <- gbm_grid$n.trees[i]
  shrink <- gbm_grid$shrinkage[i]
  min_obs <- gbm_grid$n.minobsinnode[i]
  
  # Store C-index for each fold
  c_index_fold <- numeric(cv_folds)
  
  # Inner loop for cross-validation 
  for (fold in 1:cv_folds) {
    # Split data for this fold
    train_data <- data_gbm[folds != fold, ]
    valid_data <- data_gbm[folds == fold, ]
    
    # Fit the GBM model 
    gbm_model <- gbm(SurvObj ~ ., 
                     data = train_data[, c(features, "SurvObj")],
                     distribution = "coxph", 
                     n.trees = trees, 
                     interaction.depth = depth, 
                     shrinkage = shrink,
                     n.minobsinnode = min_obs, 
                     verbose = FALSE)
    
    # Predict risk scores (using all trees for now)
    risk_scores <- predict(gbm_model, newdata = valid_data[, features], n.trees = trees) 
    
    # Calculate C-index for this fold
    c_index_fold[fold] <- concordance.index(x = risk_scores, 
                                            surv.time = valid_data$tpfs, 
                                            surv.event = valid_data$pfs)$c.index
  } 
  
  # Store average C-index for this set of hyperparameters
  results <- rbind(results, data.frame(interaction.depth = depth,
                                       n.trees = trees,
                                       shrinkage = shrink,
                                       n.minobsinnode = min_obs,
                                       c_index = mean(c_index_fold)))
}

# Print results (sorted by C-index)
results <- results[order(-results$c_index), ]
print(results)

# Get best hyperparameters
best_params <- results[1, ] 

# Fit final model using best hyperparameters (and include cv.folds)
final_gbm_model <- gbm(SurvObj ~ ., 
                       data = data_gbm[, c(features, "SurvObj")], 
                       distribution = "coxph", 
                       n.trees = best_params$n.trees, 
                       interaction.depth = best_params$interaction.depth, 
                       shrinkage = best_params$shrinkage,
                       n.minobsinnode = best_params$n.minobsinnode, 
                       cv.folds = 5,              
                       verbose = FALSE)

# Find the optimal number of trees for the final model

n_trees_optimal <- gbm.perf(final_gbm_model, method = "cv")

print(paste("Optimal number of trees for final model:", n_trees_optimal))

# Calculate the C-index for the final model 
risk_scores_final <- predict(final_gbm_model, newdata = data_gbm[, features], n.trees = n_trees_optimal)
c_index_final <- concordance.index(x = risk_scores_final, 
                                   surv.time = data_gbm$tpfs, 
                                   surv.event = data_gbm$pfs)

print(paste("C-index for the final GBM model:", c_index_final$c.index))


#--------------------------------------------------

# Save the final model and n_trees_optimal together
save(final_gbm_model, n_trees_optimal, file = "final_gbm_model.RData")

#--------------------------------------------------

# Later, load the model
load("final_gbm_model.RData")

#--------------------------------------------------
# Use the loaded model for predictions
risk_scores_final <- predict(final_gbm_model, newdata = data_gbm[, features], n.trees = n_trees_optimal)

c_index_final <- concordance.index(x = risk_scores_final, 
                                   surv.time = data_gbm$tpfs, 
                                   surv.event = data_gbm$pfs)

print(paste("C-index for the final GBM model:", c_index_final$c.index))


#--------------------------------------------------


######################################################################
######################################################################


# PREDICTIONS
# PREDICTIONS

# Create a new dataset with different levels of b_leader
new_data <- data.frame(
  b_leader = c(0, 1),
  # Add other predictors with specific values
  dnrage = rep(median(data_gbm$dnrage, na.rm = TRUE), 2),
  drb1_gvh = rep(1, 2),  # Assigning specific value
  dqb1_gvh = rep(1, 2),   # Assigning specific value
  tce = rep(0, 2),       # Assigning specific value
  condint = rep(1, 2),   # Assigning specific value
  disease_stage_rsm = rep(1, 2),  # Assigning specific value
  hctcigp = rep(0, 2),   # Assigning specific value
  rsex = rep(1, 2),      # Assigning specific value
  donorcmv = rep(1, 2),  # Assigning specific value
  recipientcmv = rep(1, 2),  # Assigning specific value
  ragecat_rsm = rep("30-39", 2),  # Assigning specific value
  relative_rsm = rep(1, 2)  # Assigning specific value
  # Add other predictors as needed
)

# Print the new dataset
print(new_data)

# Predict the linear predictor (log hazard) for the new data
log_hazard <- predict(final_gbm_model, newdata = new_data, n.trees = n_trees_optimal, type = "link")

# Convert log hazard to survival probabilities at 12, 24, and 36 months
time_points <- c(12, 24, 36)
survival_probabilities <- matrix(NA, nrow = length(time_points), ncol = length(log_hazard))

for (i in seq_along(time_points)) {
  survival_probabilities[i, ] <- exp(-exp(log_hazard) * time_points[i] / 12)  # Adjusting for monthly time points
}

# Combine the results into a data frame
results <- data.frame(
  b_leader = new_data$b_leader,
  `12_months` = survival_probabilities[1, ],
  `24_months` = survival_probabilities[2, ],
  `36_months` = survival_probabilities[3, ]
)

print(results)




######################################################################
######################################################################
######################################################################
######################################################################

library(gbm)
library(survival)
library(splines) # For spline regression

# Function to predict survival probabilities using risk scores
predict_survival_riskscore <- function(variable_level, variable_name, model, data, time_points) {
  # Create new data with the specified variable level
  new_data <- data 
  new_data[, variable_name] <- factor(variable_level, levels = levels(data[, variable_name]))
  
  # Predict risk scores
  risk_scores <- predict(model, newdata = new_data, n.trees = n_trees_optimal, type = "link")
  
  # Create survival object
  surv_obj <- Surv(data$tpfs, data$pfs)
  
  # Fit a Kaplan-Meier curve stratified by risk group
  risk_group <- cut(risk_scores, breaks = quantile(risk_scores, probs = seq(0, 1, 0.2)), include.lowest = TRUE)
  surv_fit <- survfit(surv_obj ~ risk_group, data = data.frame(surv_obj, risk_group))
  
  # Get survival probabilities at specified time points ONLY
  surv_probs_summary <- summary(surv_fit, times = time_points) 
  surv_probs <- surv_probs_summary$surv[match(time_points, surv_probs_summary$time)] 
  
  # Check if surv_probs is a matrix (multiple risk groups)
  if (is.matrix(surv_probs)) {
    # Average survival probabilities across risk groups
    return(colMeans(surv_probs))
  } else {
    # If only one risk group, return the survival probabilities directly
    return(surv_probs)
  }
}

# Time points of interest
time_points <- c(12, 24, 36)

# Loop through all variables in the 'features' vector
for (variable_name in features) {
  # Get the levels of the current variable
  variable_levels <- levels(data_gbm[, variable_name])
  
  # Loop through each level of the variable
  for (level in variable_levels) {
    surv_probs <- predict_survival_riskscore(level, variable_name, final_gbm_model, data_gbm, time_points)
    cat("Survival probabilities for", variable_name, "=", level, "at 12, 24, and 36 months:", 
        round(surv_probs, 3), "\n")
  }
  cat("\n") # Add a blank line between variables for readability
} 

# ----------------------------------

# Function to predict survival probabilities using risk score mapping
predict_survival_mapping <- function(new_data, model, time_points, data_for_mapping) {
  # Predict risk score for the new data point
  risk_score <- predict(model, newdata = new_data, n.trees = n_trees_optimal, type = "link")
  
  # Create survival object for the data used for mapping (data_gbm)
  surv_obj <- Surv(data_for_mapping$tpfs, data_for_mapping$pfs)
  
  # Fit a spline regression model for each time point
  spline_models <- lapply(time_points, function(t) {
    surv_at_t <- as.numeric(surv_obj[, "time"] >= t) # Event indicator at time t
    lm(surv_at_t ~ ns(risk_scores_original_data, df = 4), data = data_for_mapping)
  })
  
  # Predict survival probabilities for the new risk score using the fitted models
  surv_probs <- sapply(spline_models, function(m) predict(m, newdata = data.frame(risk_scores_original_data = risk_score)))
  return(surv_probs)
}

# Time points of interest
time_points <- c(12, 24, 36)

# Dnrage values to predict for
dnrage_values <- c(65, 50, 35, 20)

# Create a new data frame with representative values 
# Make sure to fill in meaningful values!
new_data_for_prediction <- data.frame(
  ragecat_rsm = factor("<=18", levels = levels(data_gbm$ragecat_rsm)), 
  rsex = factor("Male", levels = levels(data_gbm$rsex)),              
  donorcmv = factor(0, levels = levels(data_gbm$donorcmv)),           
  recipientcmv = factor(0, levels = levels(data_gbm$recipientcmv)),     
  relative_rsm = factor("1", levels = levels(data_gbm$relative_rsm)),  
  hctcigp = factor("0", levels = levels(data_gbm$hctcigp)),          
  b_leader = factor("0", levels = levels(data_gbm$b_leader)),        
  drb1_gvh = factor("0", levels = levels(data_gbm$drb1_gvh)),         
  dqb1_gvh = factor("0", levels = levels(data_gbm$dqb1_gvh)),         
  tce = factor("0", levels = levels(data_gbm$tce)),  
  disease_stage_rsm = factor("Stage1", levels = levels(data_gbm$disease_stage_rsm))  # Include disease_stage_rsm
)

# Predict risk scores for the original data (for spline model fitting)
risk_scores_original_data <- predict(final_gbm_model, newdata = data_gbm, n.trees = n_trees_optimal, type = "link")

# Loop through the dnrage values
for (dnrage_val in dnrage_values) {
  # Add dnrage to the new data frame
  new_data_for_prediction$dnrage <- dnrage_val
  
  # Predict survival probabilities using the new mapping function
  surv_probs <- predict_survival_mapping(new_data_for_prediction, final_gbm_model, time_points, data_gbm)
  
  cat("Survival probabilities for dnrage =", dnrage_val, "at 12, 24, and 36 months:", 
      round(surv_probs, 3), "\n")
}


######################################################################
######################################################################
######################################################################


#--------------------------------------------------------------------------

