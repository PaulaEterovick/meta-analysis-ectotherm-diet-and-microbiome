# FLORA OTU-Metadata Correlation Analysis
# Script to analyze correlations between OTUs and metadata variables
# with taxonomic annotation

library("readxl")
library("dplyr")      
library("tidyr")       
library("corrplot")       
library("Hmisc") 
library("ggplot2")    
library("pheatmap")      
library("RColorBrewer") 
library("reshape2")       
library("writexl") 
library("tibble")
library("magrittr")  # For %>% operator in case it's missing


data_FLORA <- read_excel("/Users/paulaeterovick/Dokumente/ProjektDFG2024/metanalysis/manuscript_meta/merged_phyloseq_complete.xlsx", sheet = "dataFLORA")
otu_mat<- read_excel("/Users/paulaeterovick/Dokumente/ProjektDFG2024/metanalysis/manuscript_meta/merged_phyloseq_complete.xlsx", sheet = "otu_counts")
tax_mat<- read_excel("/Users/paulaeterovick/Dokumente/ProjektDFG2024/metanalysis/manuscript_meta/merged_phyloseq_complete.xlsx", sheet = "taxonomy")

#testing correlations among food variables


plot(data_FLORA$crude_protein_percent, data_FLORA$vertebrate_protein_sources)
cor.test(data_FLORA$crude_protein_percent, data_FLORA$vertebrate_protein_sources, method=c("spearman"))
#S = 225756194, p-value < 2.2e-16
#  rho = 0.324468 

plot(data_FLORA$crude_protein_percent, data_FLORA$plant_protein_sources)
cor.test(data_FLORA$crude_protein_percent, data_FLORA$plant_protein_sources, method=c("spearman"))
#S = 361482376, p-value = 0.5322
#  rho = -0.01742639  

plot(data_FLORA$crude_protein_percent, data_FLORA$insect_protein_sources)
cor.test(data_FLORA$crude_protein_percent, data_FLORA$insect_protein_sources, method=c("spearman"))
#S = 415874953, p-value = 4.525e-05
#  rho =  -0.112488

plot(data_FLORA$crude_protein_percent, data_FLORA$crude_lipid_percent)
cor.test(data_FLORA$crude_protein_percent, data_FLORA$crude_lipid_percent, method=c("spearman"))
#S = 547521823, p-value < 2.2e-16
#  rho =  -0.2916801

plot(data_FLORA$crude_protein_percent, data_FLORA$vertebrate_fat)
cor.test(data_FLORA$crude_protein_percent, data_FLORA$vertebrate_fat, method=c("spearman"))
# S = 365455720, p-value = 5.219e-07
#  rho =  -0.141755 

plot(data_FLORA$crude_protein_percent, data_FLORA$plant_fat)
cor.test(data_FLORA$crude_protein_percent, data_FLORA$plant_fat, method=c("spearman"))
#S = 117325567, p-value = 8.42e-05
#  rho =  -0.1342203 

plot(data_FLORA$crude_protein_percent, data_FLORA$insect_fat)
cor.test(data_FLORA$crude_protein_percent, data_FLORA$insect_fat, method=c("spearman"))
#S = 212695264, p-value = 0.9503
#  rho =  -0.001893579 

plot(data_FLORA$crude_protein_percent, data_FLORA$total_carbohydrate_percent)
cor.test(data_FLORA$crude_protein_percent, data_FLORA$total_carbohydrate_percent, method=c("spearman"))
#S = 37571405, p-value < 2.2e-16
#  rho =  -0.5758149 

# plot(data_FLORA$crude_protein_percent, data_FLORA$total_fiber_percent)
# cor.test(data_FLORA$crude_protein_percent, data_FLORA$total_fiber_percent, method=c("spearman"))
#S = 2251416, p-value = 0.09523
#  rho =  -0.1102779

plot(data_FLORA$vertebrate_protein_sources, data_FLORA$plant_protein_sources) #maybe inverse relationship
cor.test(data_FLORA$vertebrate_protein_sources, data_FLORA$plant_protein_sources, method=c("spearman"))
#S = 490980208, p-value < 2.2e-16
#  rho =  -0.5488223

plot(data_FLORA$vertebrate_protein_sources, data_FLORA$insect_protein_sources)
cor.test(data_FLORA$vertebrate_protein_sources, data_FLORA$insect_protein_sources, method=c("spearman"))
#S = 430993014, p-value < 2.2e-16
#  rho =  -0.2896638

plot(data_FLORA$vertebrate_protein_sources, data_FLORA$crude_lipid_percent)
cor.test(data_FLORA$vertebrate_protein_sources, data_FLORA$crude_lipid_percent, method=c("spearman"))
#S = 300661030, p-value = 0.1866
#  rho =  0.03763582

plot(data_FLORA$vertebrate_protein_sources, data_FLORA$vertebrate_fat)
cor.test(data_FLORA$vertebrate_protein_sources, data_FLORA$vertebrate_fat, method=c("spearman"))
# S = 300376587, p-value = 0.02997
#  rho =  0.06156494

plot(data_FLORA$vertebrate_protein_sources, data_FLORA$plant_fat)
cor.test(data_FLORA$vertebrate_protein_sources, data_FLORA$plant_fat, method=c("spearman"))
#S = 102795343, p-value = 0.0862
#  rho =  -0.05941444

plot(data_FLORA$vertebrate_protein_sources, data_FLORA$insect_fat)
cor.test(data_FLORA$vertebrate_protein_sources, data_FLORA$insect_fat, method=c("spearman"))
# S = 196522405, p-value = 0.05181
#  rho =  -0.06043599

plot(data_FLORA$vertebrate_protein_sources, data_FLORA$total_carbohydrate_percent)
cor.test(data_FLORA$vertebrate_protein_sources, data_FLORA$total_carbohydrate_percent, method=c("spearman"))
#S = 20948486, p-value = 0.2778
#  rho =  -0.0489745

# plot(data_FLORA$vertebrate_protein_sources, data_FLORA$total_fiber_percent)
# cor.test(data_FLORA$vertebrate_protein_sources, data_FLORA$total_fiber_percent, method=c("spearman"))

#S = 623646, p-value = 7.464e-09
#  rho =  0.4090007

plot(data_FLORA$plant_protein_sources, data_FLORA$insect_protein_sources)
cor.test(data_FLORA$plant_protein_sources, data_FLORA$insect_protein_sources, method=c("spearman"))
#S = 428521257, p-value = 8.178e-14
#  rho =  -0.2061137 

plot(data_FLORA$plant_protein_sources, data_FLORA$crude_lipid_percent)
cor.test(data_FLORA$plant_protein_sources, data_FLORA$crude_lipid_percent, method=c("spearman"))
#S = 361208967, p-value = 0.002255
#  rho =  -0.08600758 

plot(data_FLORA$plant_protein_sources, data_FLORA$vertebrate_fat)
cor.test(data_FLORA$plant_protein_sources, data_FLORA$vertebrate_fat, method=c("spearman"))
#S = 319088125, p-value = 0.07063
#  rho =  -0.05175647 

plot(data_FLORA$plant_protein_sources, data_FLORA$plant_fat)
cor.test(data_FLORA$plant_protein_sources, data_FLORA$vertebrate_fat, method=c("spearman"))
# S = 319088125, p-value = 0.07063
#  rho =  -0.05175647

plot(data_FLORA$plant_protein_sources, data_FLORA$insect_fat)
cor.test(data_FLORA$plant_protein_sources, data_FLORA$insect_fat, method=c("spearman"))
#S = 204543046, p-value = 0.4228
#  rho =  -0.02461994

plot(data_FLORA$plant_protein_sources, data_FLORA$total_carbohydrate_percent)
cor.test(data_FLORA$plant_protein_sources, data_FLORA$total_carbohydrate_percent, method=c("spearman"))
#S = 20027843, p-value = 0.0002389
#  rho =  0.1599949 

# plot(data_FLORA$plant_protein_sources, data_FLORA$total_fiber_percent)
# cor.test(data_FLORA$plant_protein_sources, data_FLORA$total_fiber_percent, method=c("spearman"))
#S = 1873793, p-value = 1.117e-14
#  rho =  -0.5162836

plot(data_FLORA$insect_protein_sources, data_FLORA$crude_lipid_percent)
cor.test(data_FLORA$insect_protein_sources, data_FLORA$crude_lipid_percent, method=c("spearman"))
#S = 313019415, p-value = 0.0001333
#  rho =  0.1065394 

plot(data_FLORA$insect_protein_sources, data_FLORA$vertebrate_fat)
cor.test(data_FLORA$insect_protein_sources, data_FLORA$vertebrate_fat, method=c("spearman"))
#S = 310844293, p-value = 0.3093
#  rho =  0.02886179 

plot(data_FLORA$insect_protein_sources, data_FLORA$plant_fat)
cor.test(data_FLORA$insect_protein_sources, data_FLORA$plant_fat, method=c("spearman"))
#S = 112777324, p-value = 0.008354
#  rho =  -0.09025111 

plot(data_FLORA$insect_protein_sources, data_FLORA$insect_fat)
cor.test(data_FLORA$insect_protein_sources, data_FLORA$insect_fat, method=c("spearman"))
#S = 216535464, p-value = 0.511
#  rho =  -0.01998271 

plot(data_FLORA$insect_protein_sources, data_FLORA$total_carbohydrate_percent)
cor.test(data_FLORA$insect_protein_sources, data_FLORA$total_carbohydrate_percent, method=c("spearman"))
#S = 24020628, p-value = 0.8647
#  rho =  -0.007470025 

# plot(data_FLORA$insect_protein_sources, data_FLORA$total_fiber_percent)
# cor.test(data_FLORA$insect_protein_sources, data_FLORA$total_fiber_percent, method=c("spearman"))
#S = 1839782, p-value = 0.2388
#  rho =  -0.08030783

plot(data_FLORA$crude_lipid_percent, data_FLORA$vertebrate_fat) #maybe correlated 
cor.test(data_FLORA$crude_lipid_percent, data_FLORA$vertebrate_fat, method=c("spearman"))
#S = 117979082, p-value < 2.2e-16
#  rho =  0.6053358

plot(data_FLORA$crude_lipid_percent, data_FLORA$plant_fat) #maybe correlated
cor.test(data_FLORA$crude_lipid_percent, data_FLORA$plant_fat, method=c("spearman"))
#S = 43630262, p-value < 2.2e-16
#  rho =  0.533794 

plot(data_FLORA$crude_lipid_percent, data_FLORA$insect_fat)
cor.test(data_FLORA$crude_lipid_percent, data_FLORA$insect_fat, method=c("spearman"))
#S = 208780673, p-value = 0.8556
#  rho =  -0.005553319 

plot(data_FLORA$crude_lipid_percent, data_FLORA$total_carbohydrate_percent)
cor.test(data_FLORA$crude_lipid_percent, data_FLORA$total_carbohydrate_percent, method=c("spearman"))
#S = 10386994, p-value < 2.2e-16
#  rho =  0.56435 

# plot(data_FLORA$crude_lipid_percent, data_FLORA$total_fiber_percent)
# cor.test(data_FLORA$crude_lipid_percent, data_FLORA$total_fiber_percent, method=c("spearman"))
# S = 1613059, p-value = 0.001821
#  rho =  0.2045254

plot(data_FLORA$vertebrate_fat, data_FLORA$plant_fat) #maybe inversely correlated
cor.test(data_FLORA$vertebrate_fat, data_FLORA$plant_fat, method=c("spearman"))
#S = 95523068, p-value = 0.1454
#  rho =  -0.05097889

plot(data_FLORA$vertebrate_fat, data_FLORA$insect_fat)
cor.test(data_FLORA$vertebrate_fat, data_FLORA$insect_fat, method=c("spearman"))
#S = 167738700, p-value = 0.1423
#  rho =  0.04601486

plot(data_FLORA$vertebrate_fat, data_FLORA$total_carbohydrate_percent)
cor.test(data_FLORA$vertebrate_fat, data_FLORA$total_carbohydrate_percent, method=c("spearman"))
#S = 10036184, p-value < 2.2e-16
#  rho =  0.4974481 

# plot(data_FLORA$vertebrate_fat, data_FLORA$total_fiber_percent)
# cor.test(data_FLORA$vertebrate_fat, data_FLORA$total_fiber_percent, method=c("spearman"))
#S = 818156, p-value = 0.002108
#  rho =  0.2246733

plot(data_FLORA$plant_fat, data_FLORA$insect_fat)
cor.test(data_FLORA$plant_fat, data_FLORA$insect_fat, method=c("spearman"))
#S = NA, p-value = NA
#  rho =  NA

plot(data_FLORA$plant_fat, data_FLORA$total_carbohydrate_percent)
cor.test(data_FLORA$plant_fat, data_FLORA$total_carbohydrate_percent, method=c("spearman"))
#S = 7313279, p-value = 0.548
#  rho =  -0.03226377

# plot(data_FLORA$plant_fat, data_FLORA$total_fiber_percent)
# cor.test(data_FLORA$plant_fat, data_FLORA$total_fiber_percent, method=c("spearman"))
#S = 999671, p-value = 0.1856
#  rho =  -0.100232

plot(data_FLORA$insect_fat, data_FLORA$total_carbohydrate_percent)
cor.test(data_FLORA$insect_fat, data_FLORA$total_carbohydrate_percent, method=c("spearman"))
#S = 10062768, p-value = 0.9167
#  rho =  0.005295943

# plot(data_FLORA$insect_fat, data_FLORA$total_fiber_percent)
# cor.test(data_FLORA$insect_fat, data_FLORA$total_fiber_percent, method=c("spearman"))
#S = NA, p-value = NA
#  rho =  NA

# plot(data_FLORA$total_carbohydrate_percent, data_FLORA$total_fiber_percent)
# cor.test(data_FLORA$total_carbohydrate_percent, data_FLORA$total_fiber_percent, method=c("spearman"))
#S = 22532, p-value = 1.69e-11
#  rho =  -0.8257732




# Assume first column of otu_counts contains OTU IDs
otu_ids <- otu_mat[[1]]
otu_data <- otu_mat[, -1]

# Convert OTU data to numeric matrix
otu_matrix <- as.matrix(otu_data)
rownames(otu_matrix) <- otu_ids

# Assume first column of metadata contains sample IDs
sample_ids <- data_FLORA[[1]]
metadata_vars <- data_FLORA[, -1]

# Extract numeric metadata variables only (keep sample IDs as rownames)
numeric_cols <- sapply(metadata_vars, is.numeric)
numeric_metadata <- metadata_vars[numeric_cols]
rownames(numeric_metadata) <- sample_ids  # Explicitly set rownames to sample IDs

cat("Numeric metadata variables found:", ncol(numeric_metadata), "\n")
cat("Numeric metadata columns:", colnames(numeric_metadata), "\n")

common_samples <- intersect(colnames(otu_matrix), rownames(numeric_metadata))
cat("Common samples between OTU and metadata:", length(common_samples), "\n")

# Subset data to common samples
otu_matrix_filtered <- otu_matrix[, common_samples]
metadata_filtered <- numeric_metadata[common_samples, ]

# Remove OTUs with zero variance or very low abundance
otu_sums <- rowSums(otu_matrix_filtered)
otu_vars <- apply(otu_matrix_filtered, 1, var)

# Filter OTUs: remove those with sum < 10 or variance = 0
min_abundance <- 10
otu_keep <- otu_sums >= min_abundance & otu_vars > 0

otu_matrix_final <- otu_matrix_filtered[otu_keep, ]
otu_ids_final <- otu_ids[otu_keep]

cat("OTUs after filtering:", nrow(otu_matrix_final), "(removed", sum(!otu_keep), "low abundance/zero variance OTUs)\n")

# Transpose OTU matrix for correlation analysis (samples as rows, OTUs as columns)
otu_transposed <- t(otu_matrix_final)

# Initialize results data frame
correlation_results <- data.frame()

# Calculate correlations for each metadata variable with each OTU
for(meta_var in colnames(metadata_filtered)) {
  cat("Processing metadata variable:", meta_var, "\n")
  
  # Get metadata values for common samples
  meta_values <- metadata_filtered[[meta_var]]
  
  # Remove samples with NA values in this metadata variable
  complete_samples <- !is.na(meta_values)
  
  if(sum(complete_samples) < 3) {
    cat("Skipping", meta_var, "- insufficient non-NA samples\n")
    next
  }
  
  meta_clean <- meta_values[complete_samples]
  otu_clean <- otu_transposed[complete_samples, ]
  
  # Calculate correlations
  correlations <- cor(meta_clean, otu_clean, method = "spearman", use = "complete.obs")  
  
  # Calculate correlations and p-values for each OTU
  temp_correlations <- numeric(ncol(otu_clean))
  temp_pvalues <- numeric(ncol(otu_clean))
  
  for(j in 1:ncol(otu_clean)) {
    otu_col <- otu_clean[, j]
    if(var(otu_col) > 0) {
      tryCatch({
        test_result <- cor.test(meta_clean, otu_col, method = "spearman", exact = FALSE)
        temp_correlations[j] <- as.numeric(test_result$estimate)
        temp_pvalues[j] <- test_result$p.value
      }, error = function(e) {
        temp_correlations[j] <- 0
        temp_pvalues[j] <- 1
      })
    } else {
      temp_correlations[j] <- 0
      temp_pvalues[j] <- 1
    }
  }
  
  # Combine results
  temp_results <- data.frame(
    metadata_variable = meta_var,
    otu_id = colnames(otu_clean),
    correlation = temp_correlations,
    p_value = temp_pvalues,
    n_samples = sum(complete_samples),
    stringsAsFactors = FALSE
  )  
  
  correlation_results <- rbind(correlation_results, temp_results)
}

# Adjust p-values for multiple testing
correlation_results$p_adjusted <- p.adjust(correlation_results$p_value, method = "fdr")

# Add significance flags
correlation_results$significant <- correlation_results$p_adjusted < 0.05
correlation_results$highly_significant <- correlation_results$p_adjusted < 0.01

cat("Total correlations calculated:", nrow(correlation_results), "\n")
cat("Significant correlations (FDR < 0.05):", sum(correlation_results$significant, na.rm = TRUE), "\n")
cat("Highly significant correlations (FDR < 0.01):", sum(correlation_results$highly_significant, na.rm = TRUE), "\n")

# Add taxonomic information
cat("\nAdding taxonomic information...\n")

# Assume first column of taxonomy contains OTU IDs
taxonomy_otu_ids <- tax_mat[[1]]
taxonomy_info <- tax_mat[, -1]

# Create a lookup table for taxonomy
taxonomy_lookup <- taxonomy_info
rownames(taxonomy_lookup) <- taxonomy_otu_ids

# Match OTU IDs between correlation results and taxonomy
correlation_results_with_taxonomy <- correlation_results
correlation_results_with_taxonomy[, colnames(taxonomy_info)] <- NA

# Add taxonomic information to correlation results
for(i in 1:nrow(correlation_results)) {
  otu_id <- correlation_results$otu_id[i]
  if(otu_id %in% rownames(taxonomy_lookup)) {
    correlation_results_with_taxonomy[i, colnames(taxonomy_info)] <- taxonomy_lookup[otu_id, ]
  }
}

# Sort results by absolute correlation coefficient
correlation_results_final <- correlation_results_with_taxonomy[
  order(abs(correlation_results_with_taxonomy$correlation), decreasing = TRUE), ]

# Create summary of significant results
significant_results <- correlation_results_final[
  correlation_results_final$significant & !is.na(correlation_results_final$significant), ]


# Create output directory if it doesn't exist
output_dir <- "/Users/paulaeterovick/Dokumente/ProjektDFG2024/metanalysis/FLORALresults"
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Save complete results
write.csv(correlation_results_final, 
          file.path(output_dir, "complete_correlation_results.csv"), 
          row.names = FALSE)

# Save significant results only
write.csv(significant_results, 
          file.path(output_dir, "significant_correlations.csv"), 
          row.names = FALSE)


# Create summary statistics
summary_stats <- correlation_results_final %>%
  group_by(metadata_variable) %>%
  summarise(
    total_otus = n(),
    significant_correlations = sum(significant, na.rm = TRUE),
    highly_significant = sum(highly_significant, na.rm = TRUE),
    max_abs_correlation = ifelse(all(is.na(correlation)), NA, max(abs(correlation), na.rm = TRUE)),
    mean_abs_correlation = mean(abs(correlation), na.rm = TRUE),
    .groups = 'drop'
  )

write.csv(summary_stats, 
          file.path(output_dir, "correlation_summary.csv"), 
          row.names = FALSE)


# 1. Correlation heatmap for top significant results
if(nrow(significant_results) > 0) {
  cat("Creating heatmap with", nrow(significant_results), "significant results\n")
  
  # Create a more robust heatmap approach
  # Get top metadata variables and their top correlations
  # Filter to only include correlations with absolute value >= 0.4
  top_correlations_per_var <- significant_results %>%
    filter(abs(correlation) >= 0.4) %>%
    group_by(metadata_variable) %>%
    slice_max(abs(correlation), n = 15) %>%  # Top 15 by absolute correlation
    ungroup()
  
  cat("Correlations >= 0.4:", sum(abs(significant_results$correlation) >= 0.4), "\n")
  cat("Selected", nrow(top_correlations_per_var), "top correlations for heatmap\n")
  
  # Ensure we have enough data
  if(nrow(top_correlations_per_var) >= 2 && 
     length(unique(top_correlations_per_var$metadata_variable)) >= 1) {
    
    # Create matrix - handle duplicates by taking the first occurrence
    heatmap_matrix <- top_correlations_per_var %>%
      distinct(otu_id, metadata_variable, .keep_all = TRUE) %>%
      select(otu_id, metadata_variable, correlation) %>%
      pivot_wider(names_from = metadata_variable, 
                  values_from = correlation, 
                  values_fill = NA) %>%
      column_to_rownames("otu_id")
    
    # Convert to numeric matrix and handle any remaining issues
    heatmap_matrix <- as.matrix(heatmap_matrix)
    heatmap_matrix[is.infinite(heatmap_matrix)] <- NA
    
    # Remove rows/cols that are all NA
    heatmap_matrix <- heatmap_matrix[rowSums(!is.na(heatmap_matrix)) > 0, , drop = FALSE]
    heatmap_matrix <- heatmap_matrix[, colSums(!is.na(heatmap_matrix)) > 0, drop = FALSE]
    
    cat("Final heatmap matrix dimensions:", nrow(heatmap_matrix), "x", ncol(heatmap_matrix), "\n")
    
    # Only proceed if we have a reasonable matrix size
    if(nrow(heatmap_matrix) >= 2 && ncol(heatmap_matrix) >= 1) {
      
      # Create PDF with proper error handling - adjust dimensions based on data size
      pdf_width <- max(8, min(20, 2 + ncol(heatmap_matrix) * 0.6))
      pdf_height <- max(6, min(30, 2 + nrow(heatmap_matrix) * 0.15))
      cat("PDF dimensions:", pdf_width, "x", pdf_height, "\n")
      
      # Skip pheatmap and use ggplot2 directly for more reliable output
      cat("Creating ggplot2 heatmap for better compatibility...\n")
      
      # Create ggplot2 heatmap directly
      heatmap_success <- tryCatch({
        # Calculate dimensions - make narrower by reducing width per column
        alt_pdf_width <- max(5, min(10, 2 + length(unique(top_correlations_per_var$metadata_variable)) * 1.2))
        alt_pdf_height <- max(6, min(30, 2 + length(unique(top_correlations_per_var$otu_id)) * 0.15))
        cat("PDF dimensions:", alt_pdf_width, "x", alt_pdf_height, "\n")
        
        pdf(file.path(output_dir, "correlation_heatmap.pdf"), width = alt_pdf_width, height = alt_pdf_height)
        
        # Convert matrix back to long format for ggplot
        heatmap_long <- heatmap_matrix %>%
          as.data.frame() %>%
          rownames_to_column("otu_id") %>%
          pivot_longer(cols = -otu_id, names_to = "metadata_variable", values_to = "correlation") %>%
          filter(!is.na(correlation))
        
        # Order OTUs by metadata variable for better visualization
        heatmap_long <- heatmap_long %>%
          arrange(metadata_variable, desc(abs(correlation))) %>%
          mutate(otu_id = factor(otu_id, levels = unique(otu_id)))
        
        p <- ggplot(heatmap_long, aes(x = metadata_variable, y = otu_id, fill = correlation)) +
          geom_tile(color = "white") +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                              midpoint = 0, na.value = "grey90") +
          scale_y_discrete(position = "right") +  # Move y-axis labels to the right
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Larger x-axis text
                axis.text.y = element_text(size = 10),  # Larger y-axis text
                axis.title.x = element_text(size = 14, face = "bold"),  # Larger x-axis title
                axis.title.y = element_text(size = 14, face = "bold"),  # Larger y-axis title
                plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Larger plot title
                legend.title = element_text(size = 12, face = "bold"),  # Larger legend title
                legend.text = element_text(size = 10)) +  # Larger legend text
          labs(title = "OTU-Metadata Correlations Heatmap (|r| >= 0.4)",
               x = "Metadata Variable",
               y = "OTU ID",
               fill = "Correlation")
        
        print(p)
        dev.off()
        cat("ggplot2 heatmap saved successfully!\n")
        TRUE
        
      }, error = function(e) {
        dev.off()  # Make sure to close the PDF device
        cat("Error creating heatmap:", conditionMessage(e), "\n")
        FALSE
      })
      
    } else {
      cat("Matrix too small for heatmap (dimensions:", nrow(heatmap_matrix), "x", ncol(heatmap_matrix), ")\n")
      cat("No correlations >= 0.4 found or matrix too small. Creating dot plot instead...\n")
    }
  } else {
    cat("Insufficient data for heatmap (need at least 2 rows, found", nrow(top_correlations_per_var), ")\n")
    cat("This may be because no correlations >= 0.4 were found. Creating dot plot...\n")
  }
  
  # Always create a dot plot as backup/alternative
  if(!file.exists(file.path(output_dir, "correlation_heatmap.pdf")) || 
     nrow(top_correlations_per_var) < 4) {
    
    # Adjust dot plot dimensions based on data
    dot_width <- max(8, min(18, 3 + length(unique(significant_results$metadata_variable)) * 2))
    dot_height <- max(6, min(24, 2 + min(50, nrow(significant_results)) * 0.15))
    
    pdf(file.path(output_dir, "correlation_dotplot.pdf"), width = dot_width, height = dot_height)
    
    top_for_plot <- significant_results %>%
      slice_max(abs(correlation), n = 50)
    
    p <- ggplot(top_for_plot, aes(x = metadata_variable, y = reorder(otu_id, correlation), 
                                  color = correlation, size = abs(correlation))) +
      geom_point(alpha = 0.8) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      scale_size_continuous(range = c(1, 4)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 6)) +
      labs(title = "Top 50 OTU-Metadata Correlations",
           x = "Metadata Variable",
           y = "OTU ID",
           color = "Correlation",
           size = "Absolute Correlation")
    
    print(p)
    dev.off()
    cat("Correlation dot plot saved!\n")
  }
  
} else {
  cat("No significant results found for heatmap generation\n")
}

# 2. Distribution of correlation coefficients
pdf(file.path(output_dir, "correlation_distribution.pdf"), width = 10, height = 6)
ggplot(correlation_results_final, aes(x = correlation)) +
  geom_histogram(bins = 50, fill = "lightblue", alpha = 0.7) +
  facet_wrap(~metadata_variable, scales = "free_y") +
  theme_minimal() +
  labs(title = "Distribution of Correlation Coefficients",
       x = "Spearman Correlation Coefficient",
       y = "Frequency")
dev.off()

# 3. P-value distribution
pdf(file.path(output_dir, "pvalue_distribution.pdf"), width = 10, height = 6)
ggplot(correlation_results_final, aes(x = p_value)) +
  geom_histogram(bins = 50, fill = "lightgreen", alpha = 0.7) +
  facet_wrap(~metadata_variable, scales = "free_y") +
  theme_minimal() +
  labs(title = "Distribution of P-values",
       x = "P-value",
       y = "Frequency")
dev.off()

# 4. Volcano plot
if(length(unique(correlation_results_final$metadata_variable)) <= 4) {
  pdf(file.path(output_dir, "volcano_plots.pdf"), width = 12, height = 8)
  correlation_results_final$neg_log_p <- -log10(correlation_results_final$p_adjusted)
  
  ggplot(correlation_results_final, aes(x = correlation, y = neg_log_p)) +
    geom_point(aes(color = significant), alpha = 0.6) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
    facet_wrap(~metadata_variable) +
    theme_minimal() +
    labs(title = "Volcano Plots: Correlation vs Significance",
         x = "Spearman Correlation Coefficient",
         y = "-log10(Adjusted P-value)",
         color = "Significant") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")
  dev.off()
}

# Print summary
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("ANALYSIS COMPLETE\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Summary:\n")
print(summary_stats)
cat("\nFiles saved in directory:", output_dir, "\n")
cat("- complete_correlation_results.csv: All correlation results with taxonomy\n")
cat("- significant_correlations.csv: Only significant correlations\n")
cat("- correlation_summary.csv: Summary statistics by metadata variable\n")
if(file.exists(file.path(output_dir, "correlation_heatmap.pdf"))) {
  cat("- correlation_heatmap.pdf: Heatmap of top correlations\n")
}
if(file.exists(file.path(output_dir, "correlation_dotplot.pdf"))) {
  cat("- correlation_dotplot.pdf: Dot plot of top correlations\n")
}
cat("- correlation_distribution.pdf: Distribution of correlation coefficients\n")
cat("- pvalue_distribution.pdf: Distribution of p-values\n")
if(length(unique(correlation_results_final$metadata_variable)) <= 4) {
  cat("- volcano_plots.pdf: Volcano plots showing correlation vs significance\n")
}

# Display top 10 significant correlations
if(nrow(significant_results) > 0) {
  cat("\nTop 10 significant correlations:\n")
  print(significant_results[1:min(10, nrow(significant_results)), 
                            c("metadata_variable", "otu_id", "correlation", "p_adjusted")])
}

