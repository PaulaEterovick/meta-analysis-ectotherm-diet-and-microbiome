# Recursive Partitioning Analysis for Shannon Diversity
# Using partykit package
# Author: Analysis script
# Date: 2025-10-21

# Load required libraries
library(partykit)
library(readxl)
library(dplyr)
library(modeltools)
library(ggparty)
library(ggplot2)

# Load the data
data_path <- "/Users/paulaeterovick/Dokumente/ProjektDFG2024/metanalysis/partykit.xlsx"
data <- read_excel(data_path)

# Display basic information about the dataset
cat("Dataset dimensions:", dim(data), "\n")
cat("Column names:\n")
print(names(data))

# Check for the target variable
if(!"shannon" %in% names(data)) {
  stop("Target variable 'shannon' not found in the dataset")
}

# Define explanatory variables
explanatory_vars <- c(
  "species", "genus", "family", "order", "habitat",
  "crude_protein_percent", "vertebrate_protein_sources", "plant_protein_sources", 
  "insect_protein_sources", "crude_lipid_percent", "vertebrate_fat", "plant_fat", 
  "insect_fat", "total_carbohydrate_percent", "total_fiber_percent", "sample", 
  "preservation", "extraction_kit", "amplified_region", "sequencing", "strand", 
  "duration_experiment_days", "life_stage", "trimming"
)

# Check which variables are present in the data
available_vars <- explanatory_vars[explanatory_vars %in% names(data)]
missing_vars <- explanatory_vars[!explanatory_vars %in% names(data)]

cat("\nAvailable explanatory variables:", length(available_vars), "\n")
cat("Missing variables:", length(missing_vars), "\n")
if(length(missing_vars) > 0) {
  cat("Missing variables are:", paste(missing_vars, collapse = ", "), "\n")
}

# Data preprocessing
# Remove rows with missing shannon values
data_clean <- data[!is.na(data$shannon), ]
cat("\nRows after removing missing shannon values:", nrow(data_clean), "\n")

# Select only available variables plus shannon
model_data <- data_clean[, c("shannon", available_vars)]

# Check for missing values in explanatory variables
missing_summary <- sapply(model_data, function(x) sum(is.na(x)))
cat("\nMissing values per variable:\n")
print(missing_summary)

# Convert character variables to factors for better handling
categorical_vars <- sapply(model_data, is.character)
model_data[categorical_vars] <- lapply(model_data[categorical_vars], as.factor)

# Display summary statistics
cat("\nSummary statistics:\n")
print(summary(model_data$shannon))

# Create the formula for the recursive partitioning tree
formula_string <- paste("shannon ~", paste(available_vars, collapse = " + "))
model_formula <- as.formula(formula_string)
cat("\nModel formula:", formula_string, "\n")

# Build the recursive partitioning tree using ctree (conditional inference tree)
cat("\nBuilding conditional inference tree...\n")
shannon_tree <- ctree(model_formula, data = model_data)

# Print tree summary
print(shannon_tree)

# Plot the tree using ggparty
cat("\nPlotting the tree with ggparty...\n")

# Create ggparty plot using autoplot without default edge labels
ggparty_plot <- ggparty(shannon_tree, terminal_space = 0.3) +
  geom_edge() +
  geom_edge_label(aes(label = ifelse(id %% 2 == 0, breaks_label, "")), 
                  size = 3, shift = 0.6, hjust = 0.5) +
  geom_edge_label(aes(label = ifelse(id %% 2 == 1, breaks_label, "")), 
                  size = 3, shift = 0.4, hjust = 0.5) +
  geom_node_info() +
  geom_node_label(line_list = list(aes(label = splitvar, size = 7, nudge_x=15),
                                   aes(label = paste0("n=", nodesize))),
                  line_gpar = list(list(size = 11),
                                   list(size = 9)),
                  ids = "inner") +
  geom_node_label(aes(label = id),
                  fontface = "bold",
                  ids = "all",
                  size = 4,
                  nudge_y = 0.02) +
  geom_node_plot(gglist = list(
    geom_boxplot(aes(x = "", y = shannon),
                 width = 0.6,
                 fill = "lightblue",
                 alpha = 0.7,
                 na.rm = TRUE),
    theme_bw(base_size = 9),
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()),
    ylab("Shannon")
  ),
  shared_axis_labels = TRUE) +
  labs(title = "Recursive Partitioning Tree for Shannon Diversity") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Display the plot
print(ggparty_plot)

# Save the ggparty plot
ggsave("shannon_tree_ggparty.pdf", ggparty_plot, 
       width = 20, height = 18, units = "in")
ggsave("shannon_tree_ggparty.png", ggparty_plot, 
       width = 20, height = 18, units = "in", dpi = 300)

# Alternative: Build using mob (model-based recursive partitioning)
# This is useful if you want to fit linear models in the terminal nodes
cat("\nBuilding model-based recursive partitioning tree...\n")
# Create a simple linear model fit function
lm_fit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ...) {
  lm(y ~ 1, weights = weights, offset = offset, ...)
}
shannon_mob <- mob(shannon ~ 1 | ., data = model_data, fit = lm_fit)
print(shannon_mob)

# Plot the mob tree using ggparty
cat("\nPlotting the MOB tree with ggparty...\n")

# Create ggparty plot for MOB without default edge labels
ggparty_mob_plot <- ggparty(shannon_mob, terminal_space = 0.3) +
  geom_edge() +
  geom_edge_label(aes(label = ifelse(id %% 2 == 0, breaks_label, "")), 
                  size = 2.5, shift = 0.7, hjust = 0.5) +
  geom_edge_label(aes(label = ifelse(id %% 2 == 1, breaks_label, "")), 
                  size = 2.5, shift = 0.4, hjust = 0.5) +
  geom_node_info() +
  geom_node_label(line_list = list(aes(label = splitvar),
                                   aes(label = paste0("n=", nodesize))),
                  line_gpar = list(list(size = 11),
                                   list(size = 7)),
                  ids = "inner") +
  geom_node_label(aes(label = id),
                  fontface = "bold",
                  ids = "all",
                  nudge_y = 0.02) +
  geom_node_plot(gglist = list(
    geom_boxplot(aes(x = "", y = shannon),
                 width = 0.6,
                 fill = "lightgreen",
                 alpha = 0.7,
                 na.rm = TRUE),
    theme_bw(base_size = 9),
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()),
    ylab("Shannon")
  ),
  shared_axis_labels = TRUE) +
  labs(title = "Model-based Recursive Partitioning for Shannon Diversity") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Display the plot
print(ggparty_mob_plot)

# Save the ggparty mob plot
ggsave("shannon_mob_ggparty.pdf", ggparty_mob_plot, 
       width = 20, height = 18, units = "in")
ggsave("shannon_mob_ggparty.png", ggparty_mob_plot, 
       width = 20, height = 18, units = "in", dpi = 300)

# Variable importance (for ctree)
cat("\nVariable importance analysis...\n")
# Extract variable importance from the tree structure
tree_vars <- unique(unlist(nodeapply(shannon_tree, ids = nodeids(shannon_tree), 
                                    function(n) split_node(n)$varid)))
if(length(tree_vars) > 0) {
  var_names <- names(model_data)[tree_vars + 1]  # +1 because shannon is first column
  cat("Variables used in the tree:", paste(var_names, collapse = ", "), "\n")
} else {
  cat("No splits found in the tree (constant predictions)\n")
}

# Model diagnostics
cat("\nModel diagnostics:\n")
fitted_values <- predict(shannon_tree)
residuals <- model_data$shannon - fitted_values

# Calculate R-squared equivalent
ss_res <- sum(residuals^2, na.rm = TRUE)
ss_tot <- sum((model_data$shannon - mean(model_data$shannon, na.rm = TRUE))^2, na.rm = TRUE)
r_squared <- 1 - (ss_res / ss_tot)
cat("Pseudo R-squared:", round(r_squared, 4), "\n")

# RMSE
rmse <- sqrt(mean(residuals^2, na.rm = TRUE))
cat("RMSE:", round(rmse, 4), "\n")

# Residual plots
par(mfrow = c(2, 2))
plot(fitted_values, residuals, 
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)

hist(residuals, main = "Histogram of Residuals", xlab = "Residuals")

qqnorm(residuals, main = "Q-Q Plot of Residuals")
qqline(residuals, col = "red")

plot(model_data$shannon, fitted_values,
     main = "Observed vs Predicted",
     xlab = "Observed Shannon", ylab = "Predicted Shannon")
abline(0, 1, col = "red", lty = 2)

# Reset plotting parameters
par(mfrow = c(1, 1))

# Save results
cat("\nSaving results...\n")

# Save the tree object
saveRDS(shannon_tree, file = "shannon_tree_ctree.rds")
saveRDS(shannon_mob, file = "shannon_tree_mob.rds")

# Save predictions
predicted_mob <- predict(shannon_mob)
if(length(predicted_mob) == nrow(model_data)) {
  predictions_df <- data.frame(
    observed = model_data$shannon,
    predicted_ctree = fitted_values,
    predicted_mob = predicted_mob,
    residuals_ctree = residuals,
    residuals_mob = model_data$shannon - predicted_mob
  )
} else {
  # Handle case where mob predictions have different length
  predictions_df <- data.frame(
    observed = model_data$shannon,
    predicted_ctree = fitted_values,
    residuals_ctree = residuals
  )
  cat("Note: MOB predictions excluded due to dimension mismatch\n")
}
write.csv(predictions_df, "shannon_predictions.csv", row.names = FALSE)

cat("\nAnalysis completed successfully!\n")
cat("Tree objects saved as: shannon_tree_ctree.rds, shannon_tree_mob.rds\n")
cat("Predictions saved as: shannon_predictions.csv\n")
cat("ggparty plots saved as: shannon_tree_ggparty.pdf/png, shannon_mob_ggparty.pdf/png\n")

# Print final summary
cat("\n=== FINAL SUMMARY ===\n")
cat("Target variable: shannon\n")
cat("Number of observations:", nrow(model_data), "\n")
cat("Number of explanatory variables used:", length(available_vars), "\n")
cat("Pseudo R-squared (ctree):", round(r_squared, 4), "\n")
cat("RMSE (ctree):", round(rmse, 4), "\n")