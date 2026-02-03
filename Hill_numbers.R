#Following Alberdi and Gilbert, 2019
#Calculate Hill numbers (q=2) for alpha diversity and Faith's PD with SES

library(vegan)
library(ape)
library(iNEXT)
library(picante)
library(phangorn)

setwd("/filepath/Hill/")

# ========================================
# READ DATA
# ========================================

cat("Reading data...\n")

# Read OTU table
otu_table <- read.table("/filepath/otu_table.txt", 
                        header=TRUE, row.names=1)

cat(paste("OTU table:", nrow(otu_table), "OTUs x", ncol(otu_table), "samples\n\n"))

# Read phylogenetic tree
cat("Reading phylogenetic tree...\n")
if (file.exists("/filepath/phylogenetic_tree.newick")) {
  meta.tree <- read.tree("/filepath/phylogenetic_tree.newick")
} else {
  meta.tree <- read.tree("/filepath/phylogenetic_tree.tree")
}

cat(paste("Tree loaded:", length(meta.tree$tip.label), "tips\n\n"))

# ========================================
# ALPHA DIVERSITY - Hill numbers (q=2)
# ========================================

cat("========================================\n")
cat("ALPHA DIVERSITY - Hill numbers (q=2)\n")
cat("========================================\n\n")

# Prepare data: list of abundance vectors per sample
alpha_data <- list()
for (sample_name in colnames(otu_table)) {
  abundance_vec <- as.numeric(otu_table[, sample_name])
  names(abundance_vec) <- rownames(otu_table)
  alpha_data[[sample_name]] <- abundance_vec
}

# Remove samples with <= 1 species
valid_samples <- sapply(alpha_data, function(x) sum(x > 0) > 1)
if (sum(!valid_samples) > 0) {
  cat(paste("Removing", sum(!valid_samples), "sample(s) with <= 1 OTU\n"))
  alpha_data <- alpha_data[valid_samples]
}

cat(paste("Calculating Hill numbers for", length(alpha_data), "samples...\n"))

# Calculate Hill numbers with rarefaction/extrapolation for q=2
alpha_inext <- iNEXT(alpha_data, q = 2, datatype = "abundance", 
                     knots = 40, nboot = 50)

# Extract results - for q=2, the diversity is labeled as "Simpson diversity"
alpha_results <- alpha_inext$AsyEst[alpha_inext$AsyEst$Diversity == "Simpson diversity", ]

cat("\nAlpha diversity results:\n")
print(alpha_results[, c("Assemblage", "Observed", "Estimator", "s.e.", "LCL", "UCL")])

# Save results
write.csv(alpha_results, "alpha_diversity_hill_q2.csv", row.names=FALSE)
cat("\nSaved: alpha_diversity_hill_q2.csv\n\n")

# ========================================
# PHYLOGENETIC DIVERSITY - Faith's PD with SES
# ========================================

cat("========================================\n")
cat("PHYLOGENETIC DIVERSITY\n")
cat("========================================\n\n")

# Match OTUs with tree
otu_names <- rownames(otu_table)
tree_tips <- meta.tree$tip.label
common_otus <- intersect(otu_names, tree_tips)

cat(paste("OTUs in both table and tree:", length(common_otus), "/", length(otu_names), "\n"))

if (length(common_otus) < length(otu_names)) {
  cat(paste("Note:", length(otu_names) - length(common_otus), "OTUs excluded (not in tree)\n"))
}

# Filter data
otu_table_filtered <- otu_table[common_otus, ]
meta.tree_pruned <- keep.tip(meta.tree, common_otus)

# Root tree if necessary
if (!is.rooted(meta.tree_pruned)) {
  cat("Rooting tree at midpoint...\n")
  meta.tree_pruned <- midpoint(meta.tree_pruned)
}

# Prepare data: samples as rows, OTUs as columns
phylo_data_df <- as.data.frame(t(otu_table_filtered))

# Remove samples with <= 1 species
valid_phylo_samples <- rowSums(phylo_data_df > 0) > 1
if (sum(!valid_phylo_samples) > 0) {
  cat(paste("Removing", sum(!valid_phylo_samples), "sample(s) with <= 1 OTU\n"))
  phylo_data_df <- phylo_data_df[valid_phylo_samples, ]
}

cat(paste("\nCalculating phylogenetic diversity for", nrow(phylo_data_df), "samples...\n\n"))

# Faith's PD
cat("Calculating Faith's Phylogenetic Diversity...\n")
try_result <- try(faith_pd <- pd(phylo_data_df, meta.tree_pruned, include.root=TRUE), silent=TRUE)
if (inherits(try_result, "try-error")) {
  faith_pd <- pd(phylo_data_df, meta.tree_pruned, include.root=FALSE)
}

# Standardized Effect Size of Faith's PD
cat("Calculating SES for Faith's PD (999 randomizations)...\n")
cat("This may take a few minutes...\n")
ses_pd <- ses.pd(phylo_data_df, meta.tree_pruned, null.model = "taxa.labels", 
                 runs = 999, iterations = 1000)


# Combine results
phylo_results <- data.frame(
  Assemblage = rownames(phylo_data_df),
  Species_Richness = faith_pd$SR,
  Faith_PD = faith_pd$PD,
  SES_PD_pvalue = ses_pd$pd.obs.p
)

cat("\nPhylogenetic diversity results:\n")
print(phylo_results)

# Save results
write.csv(ses_pd, "phylogenetic_diversity_FaithPD_SES.csv", row.names=FALSE)
cat("\nSaved: phylogenetic_diversity_FaithPD_SES.csv\n\n")

# or

library(readr)
write_csv(phylo_results, "phylogenetic_diversity_FaithPD_SES.csv")
write_tsv(phylo_results, "hylogenetic_diversity_FaithPD_SES.txt")

# ========================================
# SUMMARY
# ========================================

cat("========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("Output files:\n")
cat("1. alpha_diversity_hill_q2.csv\n")
cat("   - Hill numbers (q=2) with rarefaction/extrapolation\n")
cat("   - Confidence intervals account for sampling uncertainty\n\n")

cat("2. phylogenetic_diversity_FaithPD_SES.csv\n")
cat("   - Faith's PD: Sum of phylogenetic branch lengths\n")
cat("   - SES_PD: Standardized effect size (z-score)\n")
cat("   - SES_PD_pvalue: Statistical significance\n\n")

cat("Interpretation of SES_PD:\n")
cat("  - z > 0: More phylogenetically diverse than expected\n")
cat("  - z < 0: Less diverse (phylogenetic clustering)\n")
cat("  - |z| > 1.96: Significant at p < 0.05\n\n")

cat("SES metrics control for species richness and are recommended\n")
cat("for comparisons across different methodologies.\n\n")

cat("References:\n")
cat("- Alpha diversity: Chao et al. (2014) iNEXT\n")
cat("- Phylogenetic diversity: Kembel et al. (2010) Bioinformatics 26:1463-1464\n")
cat("- SES approach: Webb et al. (2002) Annu Rev Ecol Syst 33:475-505\n")
cat("========================================\n")
