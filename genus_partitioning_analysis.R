# Recursive Partitioning Analysis for Bacterial Genera
# Using partykit package
# Date: 10.12.2025

# Load required libraries
library(partykit)
library(readxl)
library(dplyr)
library(modeltools)
library(ggparty)
library(ggplot2)

# Load the data

data <- read_excel("/Users/paulaeterovick/Dokumente/ProjektDFG2024/metanalysis/manuscript_meta/partykit_alldata.xlsx", sheet = "genera")

# Display basic information about the dataset
cat("Dataset dimensions:", dim(data), "\n")
cat("Column names:\n")
print(names(data))

# Define target variables (columns 2-29: bacterial genera)
target_variables <- names(data)[2:29]
cat("\nTarget variables (genera):", length(target_variables), "\n")
print(target_variables)

# Define explanatory variables (starting from column 30 onwards)
explanatory_vars <- c(
  "species", "genus", "family", "order", "habitat",
  "crude_protein_percent", "vertebrate_protein_sources", "plant_protein_sources", 
  "insect_protein_sources", "crude_lipid_percent", "vertebrate_fat", "plant_fat", 
  "insect_fat", "total_carbohydrate_percent", "sample", 
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

# Create output directory for results
output_dir <- "genus_partitioning_results"
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}
cat("\nResults will be saved to:", output_dir, "\n")

# Loop through each target variable (genus)
for(target_var in target_variables) {
  
  cat("\n========================================\n")
  cat("Processing:", target_var, "\n")
  cat("========================================\n")
  
  # Data preprocessing
  # Remove rows with missing values for current target variable
  data_clean <- data[!is.na(data[[target_var]]), ]
  cat("Rows after removing missing values:", nrow(data_clean), "\n")
  
  # Skip if too few observations
  if(nrow(data_clean) < 10) {
    cat("WARNING: Too few observations for", target_var, "- skipping\n")
    next
  }
  
  # Select only available variables plus current target
  model_data <- data_clean[, c(target_var, available_vars)]
  
  # Check for missing values in explanatory variables
  missing_summary <- sapply(model_data, function(x) sum(is.na(x)))
  cat("\nMissing values per variable:\n")
  print(missing_summary[missing_summary > 0])
  
  # Convert character variables to factors for better handling
  categorical_vars <- sapply(model_data, is.character)
  model_data[categorical_vars] <- lapply(model_data[categorical_vars], as.factor)
  
  # Display summary statistics
  cat("\nSummary statistics for", target_var, ":\n")
  print(summary(model_data[[target_var]]))
  
  # Create the formula for the recursive partitioning tree
  formula_string <- paste(target_var, "~", paste(available_vars, collapse = " + "))
  model_formula <- as.formula(formula_string)
  cat("\nModel formula:", formula_string, "\n")
  
  # Build the recursive partitioning tree using ctree (conditional inference tree)
  cat("\nBuilding conditional inference tree...\n")
  genus_tree <- tryCatch({
    ctree(model_formula, data = model_data)
  }, error = function(e) {
    cat("ERROR building tree for", target_var, ":", e$message, "\n")
    return(NULL)
  })
  
  # Skip if tree building failed
  if(is.null(genus_tree)) {
    cat("Skipping", target_var, "due to error\n")
    next
  }
  
  # Print tree summary
  print(genus_tree)

  # Plot the tree using ggparty
  cat("\nPlotting the tree with ggparty...\n")
  
  # Create custom edge labels that properly display numeric splits
  library(grid)
  
  # Create the ggparty object
  ggp_obj <- tryCatch({
    ggparty(genus_tree, terminal_space = 0.3)
  }, error = function(e) {
    cat("ERROR creating ggparty object for", target_var, ":", e$message, "\n")
    return(NULL)
  })
  
  if(is.null(ggp_obj)) {
    cat("Skipping plots for", target_var, "\n")
    next
  }
  
  # Fix breaks_label in the ggparty object's data
  for(i in 1:nrow(ggp_obj$data)) {
    if(!is.na(ggp_obj$data$breaks_label[i])) {
      label <- as.character(ggp_obj$data$breaks_label[i])
      # Remove "NA <= NA*" or "NA > NA*" patterns and keep only the numeric value
      # Handle both with and without spaces
      label <- gsub("NA\\s*<=\\s*NA\\*\\s*", "≤ ", label)
      label <- gsub("NA\\s*>\\s*NA\\*\\s*", "> ", label)
      # Also handle simpler patterns
      label <- gsub("^NA$", "", label)
      ggp_obj$data$breaks_label[i] <- label
    }
  }
  
  # Clean genus name for title (remove g__ prefix)
  genus_name <- gsub("g__", "", target_var)
  
  # Create ggparty plot with corrected edge labels at different heights
  # Split labels based on whether they go to left or right child (using id %% 2)
  ggparty_plot <- ggp_obj +
    geom_edge() +
    # Labels for left branches (even id after split) - positioned lower
    geom_edge_label(aes(label = ifelse(id %% 2 == 0, breaks_label, "")),
                    size = 3.5, 
                    hjust = 0.5,
                    vjust = 1.2,
                    parse = FALSE) +
    # Labels for right branches (odd id after split) - positioned higher  
    geom_edge_label(aes(label = ifelse(id %% 2 == 1, breaks_label, "")),
                    size = 3.5, 
                    hjust = 0.5,
                    vjust = -0.2,
                    parse = FALSE) +
    geom_node_info() +
    geom_node_label(line_list = list(aes(label = splitvar),
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
      geom_boxplot(aes(x = "", y = .data[[target_var]]),
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
      ylab(genus_name)
    ),
    shared_axis_labels = TRUE) +
    labs(title = paste("Recursive Partitioning Tree for", genus_name)) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  # Display the plot
  print(ggparty_plot)
  
  # Save the ggparty plot
  genus_file_name <- gsub("g__", "", target_var)
  ggsave(file.path(output_dir, paste0(genus_file_name, "_tree_ggparty.pdf")), ggparty_plot, 
         width = 20, height = 18, units = "in")
  ggsave(file.path(output_dir, paste0(genus_file_name, "_tree_ggparty.png")), ggparty_plot, 
         width = 20, height = 18, units = "in", dpi = 300)
  
  # Also create a base plot version for comparison
  cat("\nCreating base plot version...\n")
  pdf(file.path(output_dir, paste0(genus_file_name, "_tree_base.pdf")), width = 20, height = 14)
  plot(genus_tree, 
       main = paste("Recursive Partitioning Tree for", genus_name, "(Base Plot)"),
       type = "simple",
       gp = gpar(fontsize = 8),
       ip_args = list(pval = TRUE, id = TRUE),
       tp_args = list(id = TRUE))
  dev.off()
  
  png(file.path(output_dir, paste0(genus_file_name, "_tree_base.png")), width = 20, height = 14, units = "in", res = 300)
  plot(genus_tree, 
       main = paste("Recursive Partitioning Tree for", genus_name, "(Base Plot)"),
       type = "simple",
       gp = gpar(fontsize = 8),
       ip_args = list(pval = TRUE, id = TRUE),
       tp_args = list(id = TRUE))
  dev.off()
  
  cat("Plots saved for", genus_name, "\n")
}

cat("\n========================================\n")
cat("Analysis complete! All results saved to:", output_dir, "\n")
cat("========================================\n")

