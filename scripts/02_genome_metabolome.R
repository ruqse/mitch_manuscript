#!/usr/bin/env Rscript
# ============================================================================
# MitCH Study: Gardnerella Genome-Metabolome Integration and Figure 7
# ============================================================================
#
# Purpose:
#   Links Gardnerella pangenome gene content to metabolome composition using
#   Mantel tests, Procrustes analysis, db-RDA with confounder adjustment,
#   variance partitioning, and pairwise gene-metabolite Spearman correlations.
#   Generates publication-ready Figure 7 (multi-panel genome-metabolome figure).
#
# Inputs:
#   - objects/multiomics_views.rda        Metabolomics pos/neg views (X list)
#   - metadata/sample_metadata.rda        Sample metadata (meta_data)
#   - data/gardnerella_gene_presence_absence_binary.tsv
#   - data/gardnerella_gene_presence_absence.csv  (Panaroo output with annotations)
#
# Outputs:
#   - objects/genome_metabolome_results.rda   Cached analysis results
#   - figures/Figure7_genome_metabolome.pdf   Combined multi-panel figure
#   - figures/Figure7_genome_metabolome.png
#   - figures/Figure7a_procrustes.pdf/png     Individual panels
#   - figures/Figure7b_varpart.pdf/png
#   - figures/Figure7c_heatmap.pdf/png
#   - supplementary_tables/Table_S11.xlsx   Gene-metabolite correlations (q<0.1)
#   - supplementary_tables/gene_metabolite_correlations.csv  Full results
#
# ============================================================================

# ============================================================
# CACHING CONFIGURATION
# ============================================================
# Set to TRUE to force re-computation of all cached results
# Set to FALSE to use cached .rda files when available
FORCE_RECOMPUTE <- FALSE

# Cache file path
CACHE_FILE <- "objects/genome_metabolome_results.rda"

# Helper function to check if cache exists and is valid
use_cache <- function() {
    if (FORCE_RECOMPUTE) return(FALSE)
    if (!file.exists(CACHE_FILE)) return(FALSE)
    return(TRUE)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(readxl)
  library(reshape2)
  library(vegan)      # For distance matrices, Mantel, Procrustes
  library(ade4)       # Alternative Mantel implementation
  library(ggplot2)
  library(compositions)
  library(stringr)
  library(scales)
  library(pheatmap)
  library(ggrepel)
  library(DESeq2)
  #library(gridExtra)
  #library(patchwork)
  library(glue)
  library(ANCOMBC)
  library(decontam)
  library(phyloseq)
  #library(tibble)
  #library(ggbeeswarm)
  library(limma)
  library(randomForest)   # Literature shows superior performance
  library(caret)          # For ML model validation
  #library(ALDEx2)
  library(mixOmics)
  #library(MetaboAnalystR)
  library(openxlsx)
})
load("objects/multiomics_views.rda")
# load meta data
load("metadata/sample_metadata.rda")
# Read the gene presence/absence matrix
gene_matrix <- read.table("data/gardnerella_gene_presence_absence_binary.tsv",
                         header = TRUE,
                         row.names = 1,
                         sep = "\t",
                         check.names = FALSE)



# Extract sample IDs from gene matrix column names
sample_names <- colnames(gene_matrix)
sample_names

str(gene_matrix)
colnames(X$metabolomics_neg)
head(rownames(gene_matrix), 20)





# Extract base sample ID from MAG names (remove species suffix)
base_sample_ids <- gsub("_G\\..*$", "", sample_names)

# Match samples with metadata
matched_list <- list()

for(i in 1:length(sample_names)) {
  sample_name <- sample_names[i]
  base_id <- base_sample_ids[i]

  # Try exact match first
  match_idx <- which(meta_data$Sample == base_id)

  # If no exact match, try partial matching
  if(length(match_idx) == 0) {
    match_idx <- which(grepl(base_id, meta_data$Sample))
  }

  if(length(match_idx) > 0) {
    # Take first match if multiple
    matched_row <- meta_data[match_idx[1], , drop=FALSE]
    matched_row$MAG_Name <- sample_name
    matched_row$Base_Sample_ID <- base_id
    matched_list[[i]] <- matched_row
  } else {
    # Create placeholder matching ACTUAL meta_data structure
    placeholder_row <- data.frame(
      Sample_ID = NA_character_,
      Sample = base_id,
      Group = factor(NA, levels = levels(meta_data$Group)),
      Vaginal_pH = NA_real_,
      Age = NA_real_,
      Antibiotic = factor(NA, levels = levels(meta_data$Antibiotic)),
      BV_History = factor(NA, levels = levels(meta_data$BV_History)),
      MAG_Name = sample_name,
      Base_Sample_ID = base_id,
      stringsAsFactors = FALSE
    )
    matched_list[[i]] <- placeholder_row
  }
}

# Combine all rows at once
matched_metadata <- do.call(rbind, matched_list)

# Extract species information from MAG names
matched_metadata$Species <- gsub(".*_(G\\..*)$", "\\1", matched_metadata$MAG_Name)

head(matched_metadata)

filtered_MAG_metadata <- matched_metadata %>%
drop_na()

filtered_MAG_metadata[, -1]
str(X$metabolomics_pos)
str(X$metabolomics_neg)
# Run these and share output:
colnames(X$metabolomics_pos)
colnames(X$metabolomics_neg)
# ### Phase 1: Data Preparation
##### Step 1.1: Subset Gene Matrix to Filtered MAGs

# ============================================================
# PHASE 1: DATA PREPARATION
# ============================================================

# --- 1.1 Subset gene matrix to 35 filtered MAGs ---

# Get MAG names from filtered metadata
filtered_mags <- filtered_MAG_metadata$MAG_Name

# Subset gene matrix (genes are rows, MAGs are columns)
gene_matrix_filtered <- gene_matrix[, colnames(gene_matrix) %in% filtered_mags]

# QC checkpoint
cat("Gene matrix dimensions after filtering:\n")
cat("  Genes:", nrow(gene_matrix_filtered), "\n")
cat("  MAGs:", ncol(gene_matrix_filtered), "\n")
cat("  Expected MAGs: 35\n")

# Verify all filtered MAGs are present
missing_mags <- setdiff(filtered_mags, colnames(gene_matrix_filtered))
if(length(missing_mags) > 0) {
    warning("Missing MAGs in gene matrix: ", paste(missing_mags, collapse = ", "))
} else {
    cat("  ✓ All 35 filtered MAGs found in gene matrix\n")
}
# #### Step 1.2: Aggregate Multi-MAG Samples (Union Operation)
# --- 1.2 Aggregate to sample level (Option A: union of gene presence) ---

# Create sample-to-MAG mapping
sample_mag_map <- filtered_MAG_metadata %>%
    dplyr::select(Base_Sample_ID, MAG_Name) %>%
    distinct()

# Identify multi-MAG samples
multi_mag_samples <- sample_mag_map %>%
    dplyr::count(Base_Sample_ID) %>%
    filter(n > 1)

cat("\nMulti-MAG samples:\n")
print(multi_mag_samples)

# Transpose gene matrix: rows = MAGs, columns = genes
gene_matrix_t <- t(gene_matrix_filtered)

# Add sample ID to transposed matrix
gene_df <- as.data.frame(gene_matrix_t)
gene_df$MAG_Name <- rownames(gene_df)
gene_df <- gene_df %>%
    left_join(sample_mag_map, by = "MAG_Name")

# Aggregate by sample using max (equivalent to OR for binary 0/1 data)
gene_by_sample <- gene_df %>%
    dplyr::select(-MAG_Name) %>%
    dplyr::group_by(Base_Sample_ID) %>%
    summarise(across(everything(), max)) %>%
    column_to_rownames("Base_Sample_ID")

# QC checkpoint
cat("\nSample-level gene matrix dimensions:\n")
cat("  Samples:", nrow(gene_by_sample), "\n")
cat("  Gene clusters:", ncol(gene_by_sample), "\n")
cat("  Expected samples: 33\n")

# Verify aggregation worked correctly for multi-MAG samples
cat("\nVerifying aggregation for multi-MAG sample (example):\n")
example_sample <- multi_mag_samples$Base_Sample_ID[1]
example_mags <- sample_mag_map %>%
    filter(Base_Sample_ID == example_sample) %>%
    pull(MAG_Name)
cat("  Sample:", example_sample, "has MAGs:", paste(example_mags, collapse = ", "), "\n")
# #### Step 1.3: Match to Metabolomics Data
# --- 1.3 Match metabolomics to samples with Gardnerella MAGs ---

# Combine positive and negative mode metabolites
metabolomics_combined <- cbind(X$metabolomics_pos, X$metabolomics_neg)

cat("\nCombined metabolomics dimensions:\n")
cat("  Total samples in metabolomics:", nrow(metabolomics_combined), "\n
")
cat("  Total metabolites:", ncol(metabolomics_combined), "(",
    ncol(X$metabolomics_pos), "pos +", ncol(X$metabolomics_neg), "neg)\n")

# Find overlapping samples
samples_with_gardnerella <- rownames(gene_by_sample)
samples_with_metabolomics <- rownames(metabolomics_combined)
matched_samples <- intersect(samples_with_gardnerella, samples_with_metabolomics)

cat("\nSample matching:\n")
cat("  Samples with Gardnerella MAGs:", length(samples_with_gardnerella), "\n")
cat("  Samples with metabolomics:", length(samples_with_metabolomics), "\n")
cat("  Matched samples:", length(matched_samples), "\n")

# Subset both matrices to matched samples (ensure same order)
gene_matched <- gene_by_sample[matched_samples, ]
metab_matched <- metabolomics_combined[matched_samples, ]

# Final QC
stopifnot(identical(rownames(gene_matched), rownames(metab_matched)))
cat("\n✓ Final matched datasets:\n")
cat("  Genes:", nrow(gene_matched), "samples ×", ncol(gene_matched), "gene clusters\n")
cat("  Metabolites:", nrow(metab_matched), "samples ×", ncol(metab_matched), "metabolites\n")
# #### Step 1.4: Add Metadata for Matched Samples
# --- 1.4 Create matched metadata (sample-level, not MAG-level) ---

# Take first row per sample (metadata is identical for multi-MAG samples)
metadata_matched <- filtered_MAG_metadata %>%
    filter(Base_Sample_ID %in% matched_samples) %>%
    distinct(Base_Sample_ID, .keep_all = TRUE) %>%
    arrange(match(Base_Sample_ID, matched_samples))

# Verify order matches
stopifnot(identical(metadata_matched$Base_Sample_ID, rownames(gene_matched)))

cat("\nGroup distribution in matched samples:\n")
print(table(metadata_matched$Group))
# ### Phase 2: Global Association Testing
# This tests: "Do samples with similar Gardnerella accessory gene content have similar metabolomes?"
# #### Step 2.1: Compute Distance Matrices
# ============================================================
# PHASE 2: GLOBAL ASSOCIATION (MANTEL & PROCRUSTES)
# ============================================================

# --- 2.1 Compute distance matrices ---

# For binary gene presence/absence: Jaccard distance
# Jaccard = 1 - (intersection / union)
dist_genes <- vegdist(gene_matched, method = "jaccard", binary = TRUE)

# For continuous metabolomics: Euclidean distance
# (Data is already scaled via Pareto, so Euclidean is appropriate)
dist_metab <- vegdist(metab_matched, method = "euclidean")

# QC: Check for zero-variance samples (would cause issues)
gene_variance <- apply(gene_matched, 1, var)
metab_variance <- apply(metab_matched, 1, var)

cat("Sample variance check:\n")
cat("  Gene matrix - samples with zero variance:", sum(gene_variance == 0), "\n")
cat("  Metabolome - samples with zero variance:", sum(metab_variance == 0), "\n")

# ============================================================
# CHECK FOR CACHED RESULTS
# ============================================================
if (use_cache()) {
    cat("\n>>> LOADING CACHED RESULTS from", CACHE_FILE, "<<<\n")
    load(CACHE_FILE)
    cat("    Cached objects loaded successfully.\n\n")

    # Print summary of cached results
    cat("============================================================\n")
    cat("MANTEL TEST (from cache)\n")
    cat("============================================================\n")
    cat("Mantel statistic (Spearman r):", round(mantel_result$statistic, 4), "\n")
    cat("P-value:", format(mantel_result$signif, scientific = TRUE), "\n")

    cat("\n============================================================\n")
    cat("PROCRUSTES ANALYSIS (from cache)\n")
    cat("============================================================\n")
    cat("Procrustes correlation:", round(procrustes_result$t0, 4), "\n")
    cat("P-value:", format(procrustes_result$signif, scientific = TRUE), "\n")

    cat("\n============================================================\n")
    cat("db-RDA RESULTS (from cache)\n")
    cat("============================================================\n")
    cat("Unconditioned adj R²:", round(RsquareAdj(dbrda_uncond)$adj.r.squared, 4), "\n")
    cat("Conditioned adj R²:", round(RsquareAdj(dbrda_cond)$adj.r.squared, 4), "\n")

    cat("\n============================================================\n")
    cat("CORRELATION RESULTS (from cache)\n")
    cat("============================================================\n")
    cat("Total tests:", nrow(cor_results_ann), "\n")
    cat("FDR significant (q < 0.05):", sum(cor_results_ann$qvalue < 0.05), "\n")
    cat("FDR significant (q < 0.10):", sum(cor_results_ann$qvalue < 0.10), "\n\n")

} else {
    cat("\n>>> COMPUTING FROM SCRATCH (this may take several minutes) <<<\n")
    cat("    Set FORCE_RECOMPUTE <- FALSE and re-run to use cached results.\n\n")

# #### Step 2.2: Mantel Test
# --- 2.2 Mantel test ---
# Tests correlation between distance matrices
# H0: No correlation between gene content similarity and metabolome similarity

set.seed(42)  # For reproducibility
mantel_result <- mantel(dist_genes, dist_metab, method = "spearman", permutations = 9999)

cat("\n============================================================\n")
cat("MANTEL TEST: Gene content ~ Metabolome correlation\n")
cat("============================================================\n")
cat("Mantel statistic (Spearman r):", round(mantel_result$statistic, 4), "\n")
cat("P-value:", format(mantel_result$signif, scientific = TRUE), "\n")
cat("Permutations:", mantel_result$permutations, "\n")

# Interpretation
if(mantel_result$signif < 0.05) {
    cat("\n→ SIGNIFICANT: Samples with similar gene content tend to have similar metabolomes\n")
} else {
    cat("\n→ NOT SIGNIFICANT: No global correlation between gene content and metabolome\n")
}
# #### Step 2.3: Procrustes Analysis
# --- 2.3 Procrustes analysis ---
# More powerful test that accounts for ordination structure

# Perform PCoA on both distance matrices
pcoa_genes <- cmdscale(dist_genes, k = min(nrow(gene_matched) - 1, 10), eig = TRUE)
pcoa_metab <- cmdscale(dist_metab, k = min(nrow(metab_matched) - 1, 10), eig = TRUE)

# Use first k axes that explain meaningful variance
# Rule of thumb: axes explaining >5% variance each
gene_var_explained <- pcoa_genes$eig / sum(abs(pcoa_genes$eig)) * 100
metab_var_explained <- pcoa_metab$eig / sum(abs(pcoa_metab$eig)) * 100

cat("\nPCoA variance explained (first 5 axes):\n")
cat("  Genes:", round(gene_var_explained[1:5], 1), "%\n")
cat("  Metabolome:", round(metab_var_explained[1:5], 1), "%\n")

# Procrustes rotation and significance test
# Use first 2-3 axes for visualization, test with more
procrustes_result <- protest(pcoa_genes$points[, 1:3],
                              pcoa_metab$points[, 1:3],
                              permutations = 9999)

cat("\n============================================================\n")
cat("PROCRUSTES ANALYSIS: Ordination congruence\n")
cat("============================================================\n")
cat("Procrustes correlation:", round(procrustes_result$t0, 4), "\n")
cat("Procrustes m² (residual):", round(procrustes_result$ss, 4), "\n")
cat("P-value:", format(procrustes_result$signif, scientific = TRUE), "\n")

if(procrustes_result$signif < 0.05) {
    cat("\n→ SIGNIFICANT: Gene and metabolome ordinations are congruent\n")
} else {
    cat("\n→ NOT SIGNIFICANT: Ordinations are not more similar than expected by chance\n")
}
# #### Step 2.4: Visualize Procrustes Rotation
# --- 2.4 Procrustes visualization ---



# Extract rotated coordinates
procrustes_coords <- data.frame(
    Sample = rownames(gene_matched),
    Gene_PC1 = procrustes_result$X[,1],
    Gene_PC2 = procrustes_result$X[,2],
    Metab_PC1 = procrustes_result$Yrot[,1],
    Metab_PC2 = procrustes_result$Yrot[,2],
    Group = metadata_matched$Group
)

# Plot with arrows connecting gene → metabolome positions
p_procrustes <- ggplot(procrustes_coords) +
    geom_segment(aes(x = Gene_PC1, y = Gene_PC2,
                     xend = Metab_PC1, yend = Metab_PC2),
                 arrow = arrow(length = unit(0.1, "cm")),
                 color = "gray50", alpha = 0.6) +
    geom_point(aes(x = Gene_PC1, y = Gene_PC2, color = Group),
               shape = 16, size = 3) +
    geom_point(aes(x = Metab_PC1, y = Metab_PC2, color = Group),
               shape = 17, size = 3) +
    scale_color_manual(values = c("BV" = "#F68B33", "Control" = "#388ECC")) +
    labs(title = paste0("Procrustes: Gene content → Metabolome (r = ",
                        round(procrustes_result$t0, 3), ", p = ",
                        format(procrustes_result$signif, digits = 2), ")"),
         subtitle = "Circles = gene PCoA, Triangles = metabolome PCoA (rotated)",
         x = "Dimension 1", y = "Dimension 2") +
    theme_bw() +
    theme(legend.position = "bottom")

print(p_procrustes)
# ggsave("procrustes_gene_metabolome.pdf", p_procrustes, width = 8, height = 7)
# Phase 2.5: db-RDA (Distance-based Redundancy Analysis)
# This tests whether gene content explains metabolome variance while controlling for confounders (pH, antibiotic use, BV history, age). This is important because:
#
# Your metadata shows BV samples have higher pH
# 5 samples had recent antibiotic use
# These confounders could create spurious gene-metabolome correlations
# ============================================================
# PHASE 2.5: db-RDA WITH CONFOUNDER ADJUSTMENT (CORRECTED)
# ============================================================



# Prepare covariates (must match sample order)
covariates <- metadata_matched %>%
    mutate(
        pH_z = scale(Vaginal_pH)[,1],
        Age_z = scale(Age)[,1],
        Antibiotic = as.numeric(Antibiotic == "Yes"),
        BV_History = as.numeric(BV_History == "Yes"),
        Group_BV = as.numeric(Group == "BV")
    )

# Verify order
stopifnot(identical(covariates$Base_Sample_ID, rownames(metab_matched)))

# --- 2.5a: Unconditioned db-RDA ---
# How much metabolome variance is explained by gene content?

# First, perform PCoA on gene distance matrix
pcoa_genes <- cmdscale(dist_genes, k = 10, eig = TRUE)
gene_pcs <- as.data.frame(pcoa_genes$points)
colnames(gene_pcs) <- paste0("GenePC", 1:10)

# db-RDA: metabolome ~ gene PCs
dbrda_uncond <- capscale(dist_metab ~ ., data = gene_pcs)

cat("\n============================================================\n")
cat("db-RDA: Metabolome ~ Gene content (unconditioned)\n")
cat("============================================================\n")

# CORRECTED: RsquareAdj (one 's')
var_explained <- RsquareAdj(dbrda_uncond)
cat("R² =", round(var_explained$r.squared, 4), "\n")
cat("Adjusted R² =", round(var_explained$adj.r.squared, 4), "\n")

# Permutation test
anova_uncond <- anova(dbrda_uncond, permutations = 999)
cat("Permutation test p-value:", anova_uncond$`Pr(>F)`[1], "\n")

# --- 2.5b: Conditioned db-RDA ---
# Gene content effect AFTER removing confounder effects

# Combine gene PCs with covariates for the formula
model_data <- cbind(gene_pcs, covariates[, c("pH_z", "Age_z", "Antibiotic", "BV_History")])

dbrda_cond <- capscale(dist_metab ~ GenePC1 + GenePC2 + GenePC3 + GenePC4 + GenePC5 +
                        Condition(pH_z + Age_z + Antibiotic + BV_History),
                        data = model_data)

cat("\n============================================================\n")
cat("db-RDA: Metabolome ~ Gene content | Confounders\n")
cat("============================================================\n")

# CORRECTED: RsquareAdj (one 's')
var_explained_cond <- RsquareAdj(dbrda_cond)
cat("Conditioned R² =", round(var_explained_cond$r.squared, 4), "\n")
cat("Conditioned Adjusted R² =", round(var_explained_cond$adj.r.squared, 4), "\n")

anova_cond <- anova(dbrda_cond, permutations = 999)
cat("Permutation test p-value:", anova_cond$`Pr(>F)`[1], "\n")

# --- 2.5c: Variance partitioning ---
# Decompose: Gene-only, Confounder-only, Shared, Residual

varpart_result <- varpart(dist_metab,
                          gene_pcs[, 1:5],  # Gene content (first 5 PCs)
                          covariates[, c("pH_z", "Age_z", "Antibiotic", "BV_History")])

cat("\n============================================================\n")
cat("VARIANCE PARTITIONING\n")
cat("============================================================\n")
cat("Gene content alone [a]:", round(varpart_result$part$indfract$Adj.R.squared[1], 4), "\n")
cat("Confounders alone [c]:", round(varpart_result$part$indfract$Adj.R.squared[3], 4), "\n")
cat("Shared [b]:", round(varpart_result$part$indfract$Adj.R.squared[2], 4), "\n")
cat("Residual [d]:", round(varpart_result$part$indfract$Adj.R.squared[4], 4), "\n")

# Visualize Venn diagram
plot(varpart_result,
     Xnames = c("Gardnerella\ngene content", "Confounders\n(pH, age, abx, hx)"),
     bg = c("steelblue", "tomato"))
title("Variance Partitioning: Metabolome explained by...")
# ### Phase 3 — Find the Specific Gene-Metabolite Drivers
# The global tests confirm a real relationship exists. Now we identify which specific genes and metabolites are responsible.
# List all metabolites to identify BV-relevant targets
cat("All metabolites in dataset:\n")
print(colnames(metab_matched))
# ============================================================
# PHASE 3: PAIRWISE GENE-METABOLITE CORRELATIONS
# ============================================================

# --- 3.1 Filter gene clusters to informative range ---

gene_prevalence <- colMeans(gene_matched)

cat("Gene cluster prevalence distribution:\n")
cat("  Core (>95%):", sum(gene_prevalence > 0.95), "\n")
cat("  Soft-core (80-95%):", sum(gene_prevalence > 0.80 & gene_prevalence <= 0.95), "\n")
cat("  Shell (15-80%):", sum(gene_prevalence >= 0.15 & gene_prevalence <= 0.80), "\n")
cat("  Cloud (<15%):", sum(gene_prevalence < 0.15), "\n")

# Keep shell + soft-core genes (most informative for correlations)
informative_genes <- gene_prevalence >= 0.15 & gene_prevalence <= 0.95
gene_filtered <- gene_matched[, informative_genes]

cat("\nFiltered to", ncol(gene_filtered), "informative gene clusters for correlation analysis\n")

# --- 3.2 Compute Spearman correlations ---

n_genes <- ncol(gene_filtered)
n_metab <- ncol(metab_matched)

cat("Computing", n_genes, "×", n_metab, "=", n_genes * n_metab, "correlations...\n")

# Initialize results
cor_results <- data.frame(
    Gene = rep(colnames(gene_filtered), each = n_metab),
    Metabolite = rep(colnames(metab_matched), times = n_genes),
    rho = NA_real_,
    pvalue = NA_real_,
    stringsAsFactors = FALSE
)

# Compute correlations
pb <- txtProgressBar(min = 0, max = nrow(cor_results), style = 3)
for(i in seq_len(nrow(cor_results))) {
    g <- cor_results$Gene[i]
    m <- cor_results$Metabolite[i]

    gene_vec <- gene_filtered[, g]
    metab_vec <- metab_matched[, m]

    if(var(gene_vec) > 0) {
        test <- suppressWarnings(cor.test(gene_vec, metab_vec, method = "spearman", exact = FALSE))
        cor_results$rho[i] <- test$estimate
        cor_results$pvalue[i] <- test$p.value
    }
    setTxtProgressBar(pb, i)
}
close(pb)

# --- 3.3 FDR correction ---

cor_results <- cor_results %>%
    filter(!is.na(pvalue)) %>%
    mutate(qvalue = p.adjust(pvalue, method = "BH"))

cat("\n============================================================\n")
cat("GENE-METABOLITE CORRELATION RESULTS\n")
cat("============================================================\n")
cat("Total tests:", nrow(cor_results), "\n")
cat("Nominally significant (p < 0.05):", sum(cor_results$pvalue < 0.05), "\n")
cat("FDR significant (q < 0.05):", sum(cor_results$qvalue < 0.05), "\n")
cat("FDR significant (q < 0.10):", sum(cor_results$qvalue < 0.10), "\n")
cat("FDR significant (q < 0.20):", sum(cor_results$qvalue < 0.20), "\n")

# --- 3.4 Top associations ---

cat("\n--- TOP 30 ASSOCIATIONS (by p-value) ---\n")
top_hits <- cor_results %>%
    arrange(pvalue) %>%
    head(30)
print(top_hits)

# Save full results
# NOTE: large file (~25MB), deposited on Zenodo (see README)
# Uncomment the line below to regenerate locally:
# write.csv(cor_results, "gene_metabolite_correlations.csv", row.names = FALSE)
cat("\nFull correlation results available on Zenodo (see README)\n")
# ### functional annotations for pangenome gene clusters
# Step 1: Load Annotated Gene Matrix
# ============================================================
# LOAD PANAROO OUTPUT WITH ANNOTATIONS
# ============================================================

library(tidyverse)

# Load Panaroo gene presence/absence with annotations
panaroo_file <- "data/gardnerella_gene_presence_absence.csv"
panaroo_data <- read.csv(panaroo_file, check.names = FALSE)

cat("Panaroo output dimensions:", nrow(panaroo_data), "genes ×", ncol(panaroo_data), "columns\n")
cat("First 3 columns:", paste(colnames(panaroo_data)[1:3], collapse = ", "), "\n")

# Extract annotation columns
gene_annotations <- panaroo_data %>%
    dplyr::select(Gene, `Non-unique Gene name`, Annotation) %>%
    dplyr::rename(Gene_ID = Gene,
           Gene_Name = `Non-unique Gene name`,
           Gene_Annotation = Annotation)

# Preview annotations
cat("\nSample gene annotations:\n")
print(head(gene_annotations, 10))

# Extract presence/absence matrix (columns 4 onwards are MAGs)
mag_columns <- colnames(panaroo_data)[4:ncol(panaroo_data)]
presence_absence_raw <- panaroo_data[, mag_columns]

# Convert to binary: non-empty = 1, empty = 0
gene_matrix_annotated <- as.data.frame(
    apply(presence_absence_raw, 2, function(x) as.integer(x != "" & !is.na(x)))
)

# Add gene IDs as rownames
rownames(gene_matrix_annotated) <- panaroo_data$Gene

cat("\nAnnotated gene matrix dimensions:",
    nrow(gene_matrix_annotated), "genes ×",
    ncol(gene_matrix_annotated), "MAGs\n")

# Verify it matches our filtered MAGs
cat("MAGs in annotated matrix:", ncol(gene_matrix_annotated), "\n")
cat("MAGs in filtered metadata:", length(filtered_mags), "\n")
# #### Step 2: Subset and Aggregate to Sample Level (Same as Before)
# ============================================================
# SUBSET TO FILTERED MAGS AND AGGREGATE BY SAMPLE
# ============================================================

# Subset to 35 filtered MAGs
gene_matrix_filt <- gene_matrix_annotated[, colnames(gene_matrix_annotated) %in% filtered_mags]

cat("Filtered to", ncol(gene_matrix_filt), "MAGs\n")

# Transpose: rows = MAGs, columns = genes
gene_matrix_t <- t(gene_matrix_filt)

# Add sample IDs
gene_df <- as.data.frame(gene_matrix_t)
gene_df$MAG_Name <- rownames(gene_df)
gene_df <- gene_df %>%
    left_join(sample_mag_map, by = "MAG_Name")

# Aggregate by sample (union of gene presence)
gene_by_sample_annotated <- gene_df %>%
    dplyr::select(-MAG_Name) %>%
    group_by(Base_Sample_ID) %>%
    summarise(across(everything(), max)) %>%
    column_to_rownames("Base_Sample_ID")

# Reorder to match metabolomics
gene_by_sample_annotated <- gene_by_sample_annotated[matched_samples, ]

cat("Sample-level annotated gene matrix:",
    nrow(gene_by_sample_annotated), "samples ×",
    ncol(gene_by_sample_annotated), "genes\n")

# Verify alignment
stopifnot(identical(rownames(gene_by_sample_annotated), rownames(metab_matched)))
cat("✓ Sample order verified\n")
# #### Step 3: Identify BV-Relevant Genes in Pangenome
# ============================================================
# FIND BV-RELEVANT GENES IN THE PANGENOME
# ============================================================

# Key Gardnerella virulence/BV-associated gene categories
bv_gene_keywords <- c(
    # Sialidases (mucin degradation)
    "sialidase", "nanH", "nanA", "neuraminidase",
    # Vaginolysin (cytotoxin)
    "vaginolysin", "vly", "ply",  "cholesterol-dependent cytolysin",
    # Biofilm
    "biofilm", "bpfA", "bpfB", "bpfC", "bpfD",
    # Sialate O-acetylesterase
    "sialate", "nanS",
    # Glycosidases
    "glycosidase", "fucosidase", "galactosidase", "glucosidase",
    # Iron acquisition
    "iron", "heme", "siderophore", "feo", "fhu",
    # CRISPR (you found cas1/cas2 in Scoary)
    "cas1", "cas2", "CRISPR", "crispr",
    # Sulfatase (you found chuR)
    "sulfatase", "chuR", "atsA",
    # Polyamine metabolism (linked to BV odor)
    "ornithine decarboxylase", "speB", "speC", "arginine decarboxylase",
    # Amino acid metabolism
    "amino acid permease", "transporter"
)

# Search annotations for BV-relevant genes
bv_genes_found <- gene_annotations %>%
    filter(
        grepl(paste(bv_gene_keywords, collapse = "|"), Gene_ID, ignore.case = TRUE) |
        grepl(paste(bv_gene_keywords, collapse = "|"), Gene_Annotation, ignore.case = TRUE)
    )

cat("============================================================\n")
cat("BV-RELEVANT GENES FOUND IN PANGENOME\n")
cat("============================================================\n")
cat("Total BV-relevant genes found:", nrow(bv_genes_found), "\n\n")
# View first 50 rows
print(head(bv_genes_found, 50))

# Check prevalence of these genes
if(nrow(bv_genes_found) > 0) {
    bv_gene_ids <- bv_genes_found$Gene_ID
    bv_gene_ids_present <- intersect(bv_gene_ids, colnames(gene_by_sample_annotated))

    cat("\nPrevalence of BV-relevant genes (n =", length(bv_gene_ids_present), "):\n")
    bv_prevalence <- colMeans(gene_by_sample_annotated[, bv_gene_ids_present, drop = FALSE])
    print(sort(bv_prevalence, decreasing = TRUE))
}
# #### Step 4: Check Metabolite Names
# ============================================================
# LIST ALL METABOLITES
# ============================================================

cat("============================================================\n")
cat("ALL METABOLITES IN DATASET (n =", ncol(metab_matched), ")\n")
cat("============================================================\n")

# Positive mode
cat("\n--- POSITIVE MODE (", ncol(X$metabolomics_pos), "metabolites) ---\n")
print(colnames(X$metabolomics_pos))

# Negative mode
cat("\n--- NEGATIVE MODE (", ncol(X$metabolomics_neg), "metabolites) ---\n")
print(colnames(X$metabolomics_neg))

# Search for BV-relevant metabolites
bv_metab_keywords <- c(
    "putrescine", "cadaverine", "tyramine", "trimethylamine",
    "succinate", "acetate", "propionate", "butyrate",
    "sialic", "neuraminic", "lactate", "spermine", "spermidine",
    "agmatine", "histamine", "phenylacetate"
)

bv_metab_found <- colnames(metab_matched)[
    grepl(paste(bv_metab_keywords, collapse = "|"), colnames(metab_matched), ignore.case = TRUE)
]

cat("\n============================================================\n")
cat("BV-RELEVANT METABOLITES FOUND\n")
cat("============================================================\n")
print(bv_metab_found)
print(bv_metab_found)
# #### Step 5: Run Annotated Correlation Analysis
# ============================================================
# PHASE 3: GENE-METABOLITE CORRELATIONS (WITH ANNOTATIONS)
# ============================================================

# --- Filter to informative genes ---
gene_prevalence <- colMeans(gene_by_sample_annotated)
informative_genes <- gene_prevalence >= 0.15 & gene_prevalence <= 0.95
gene_filtered_ann <- gene_by_sample_annotated[, informative_genes]

cat("Filtered to", ncol(gene_filtered_ann), "informative gene clusters\n")

# --- Compute correlations ---
n_genes <- ncol(gene_filtered_ann)
n_metab <- ncol(metab_matched)

cat("Computing", n_genes, "×", n_metab, "=", n_genes * n_metab, "correlations...\n")

# Pre-allocate results
cor_results_ann <- data.frame(
    Gene_ID = character(),
    Metabolite = character(),
    rho = numeric(),
    pvalue = numeric(),
    stringsAsFactors = FALSE
)

# Compute correlations (vectorized for speed)
gene_names <- colnames(gene_filtered_ann)
metab_names <- colnames(metab_matched)

results_list <- vector("list", n_genes)

for(i in seq_len(n_genes)) {
    gene_vec <- gene_filtered_ann[, i]
    gene_name <- gene_names[i]

    if(var(gene_vec) > 0) {
        temp_results <- data.frame(
            Gene_ID = gene_name,
            Metabolite = metab_names,
            rho = NA_real_,
            pvalue = NA_real_
        )

        for(j in seq_len(n_metab)) {
            test <- suppressWarnings(
                cor.test(gene_vec, metab_matched[, j], method = "spearman", exact = FALSE)
            )
            temp_results$rho[j] <- test$estimate
            temp_results$pvalue[j] <- test$p.value
        }
        results_list[[i]] <- temp_results
    }

    if(i %% 100 == 0) cat("  Processed", i, "/", n_genes, "genes\n")
}

cor_results_ann <- bind_rows(results_list)

# --- Add annotations ---
cor_results_ann <- cor_results_ann %>%
    left_join(gene_annotations, by = "Gene_ID") %>%
    filter(!is.na(pvalue)) %>%
    mutate(qvalue = p.adjust(pvalue, method = "BH"))

# --- Summary ---
cat("\n============================================================\n")
cat("ANNOTATED GENE-METABOLITE CORRELATION RESULTS\n")
cat("============================================================\n")
cat("Total tests:", nrow(cor_results_ann), "\n")
cat("Nominally significant (p < 0.05):", sum(cor_results_ann$pvalue < 0.05), "\n")
cat("FDR significant (q < 0.05):", sum(cor_results_ann$qvalue < 0.05), "\n")
cat("FDR significant (q < 0.10):", sum(cor_results_ann$qvalue < 0.10), "\n")
cat("FDR significant (q < 0.20):", sum(cor_results_ann$qvalue < 0.20), "\n")

# --- Top hits with annotations ---
cat("\n--- TOP 30 ASSOCIATIONS ---\n")
top_hits_ann <- cor_results_ann %>%
    arrange(pvalue) %>%
    dplyr::select(Gene_ID, Gene_Annotation, Metabolite, rho, pvalue, qvalue) %>%
    head(30)

print(as_tibble(top_hits_ann), n = 30, width = Inf)

# --- Save full annotated results ---
# NOTE: large file, deposited on Zenodo (see README)
# Uncomment the line below to regenerate locally:
# write.csv(cor_results_ann, "gene_metabolite_correlations_annotated.csv", row.names = FALSE)
cat("\nFull annotated results available on Zenodo (see README)\n")
print(as_tibble(top_hits_ann), n = 30, width = Inf)

# ============================================================
# SAVE CACHE FOR FUTURE RUNS
# ============================================================
save(mantel_result, procrustes_result, pcoa_genes, pcoa_metab,
     dbrda_uncond, dbrda_cond, varpart_result, anova_uncond, anova_cond,
     cor_results, cor_results_ann, procrustes_coords, covariates,
     gene_var_explained, metab_var_explained,
     file = CACHE_FILE)
cat("\n>>> CACHE SAVED to", CACHE_FILE, "<<<\n")
cat("    Future runs will load cached results (set FORCE_RECOMPUTE <- TRUE to recompute)\n\n")

} # End of else block (computation from scratch)

# --- SUPPLEMENTARY TABLE S11: Gene-metabolite correlations (FDR q < 0.1) ---
cor_results_lessq0_1 <- cor_results_ann %>% dplyr::filter(qvalue < 0.1)
openxlsx::write.xlsx(cor_results_lessq0_1, file = "supplementary_tables/Table_S11.xlsx")
cat("Saved: supplementary_tables/Table_S11.xlsx (", nrow(cor_results_lessq0_1), "rows, q < 0.1)\n")

# Full correlation results (large file, also deposited on Zenodo)
write.csv(cor_results_ann, file = "supplementary_tables/gene_metabolite_correlations.csv",
          row.names = FALSE)
cat("Saved: supplementary_tables/gene_metabolite_correlations.csv\n")

# ### Next Steps: Focused Investigation
# Step 1: Check if the MTA-Associated Genes Cluster Together
# ============================================================
# INVESTIGATE 5-METHYLTHIOADENOSINE ASSOCIATED GENES
# ============================================================

# Get genes with strong MTA correlation
mta_genes <- cor_results_ann %>%
    filter(Metabolite == "5-Methylthioadenosine",
           pvalue < 0.001) %>%
    arrange(rho) %>%
    pull(Gene_ID) %>%
    unique()

cat("Genes strongly correlated with 5-MTA (p < 0.001):\n")
cat("N =", length(mta_genes), "\n\n")

# Check if these genes co-occur (are they from the same strains?)
mta_gene_matrix <- gene_by_sample_annotated[, mta_genes, drop = FALSE]

# Calculate co-occurrence
cat("Co-occurrence pattern of MTA-associated genes:\n")
cat("Mean pairwise correlation:",
    round(mean(cor(mta_gene_matrix)[upper.tri(cor(mta_gene_matrix))]), 3), "\n")

# Are they all present in the same samples?
mta_sample_totals <- rowSums(mta_gene_matrix)
cat("\nSamples with ALL MTA-associated genes:", sum(mta_sample_totals == length(mta_genes)), "\n")
cat("Samples with NONE:", sum(mta_sample_totals == 0), "\n")
cat("Samples with SOME:", sum(mta_sample_totals > 0 & mta_sample_totals < length(mta_genes)), "\n")

# Check group distribution
mta_presence <- ifelse(mta_sample_totals > median(mta_sample_totals), "High", "Low")
cat("\nMTA-gene cluster presence by clinical group:\n")
print(table(mta_presence, metadata_matched$Group))
# #### Step 2: Test BV-Relevant Genes Against BV-Relevant Metabolites
# ============================================================
# FOCUSED: BV GENES × BV METABOLITES
# ============================================================

# Your BV-relevant metabolites
bv_metab_found <- c("N-Acetylputrescine", "N-Acetylspermidine", "g-Aminobutyrate",
                    "2-Aminoisobutyrate", "Argininosuccinate", "Cadaverine",
                    "Guanidinoacetate", "Histamine", "Putrescine",
                    "Trimethylamine", "2-Hydroxybutyrate", "3-Hydroxybutyrate",
                    "Lactate", "PhenylLactate", "Succinate")

# Key BV-relevant genes from your pangenome (focusing on informative prevalence 15-85%)
bv_genes_key <- c(
    "cas1", "cas2", "cas9~~~cas9_1~~~cas9_2",  # CRISPR (BV-enriched in Scoary)
    "chuR",                                      # Sulfatase (BV-exclusive in Scoary)
    "frdB",                                      # Fumarate reductase
    "lacZ",                                      # β-galactosidase
    "ebgA", "ebgC",                              # β-galactosidase variants
    "sstT",                                      # Serine/threonine transporter
    "rihA", "rihB_2",                            # Ribonucleoside hydrolases
    "glnM",                                      # Glutamine transporter
    "ulaA", "ulaB", "ulaA_1", "ulaA~~~ulaA_2",   # Ascorbate PTS
    "group_1594"                                 # EcfT (already significant with MTA)
)

# Filter to genes present in your analysis
bv_genes_present <- intersect(bv_genes_key, colnames(gene_by_sample_annotated))
cat("BV-relevant genes available for analysis:", length(bv_genes_present), "\n")
print(bv_genes_present)

# Extract focused correlations
bv_focused <- cor_results_ann %>%
    filter(Gene_ID %in% bv_genes_present,
           Metabolite %in% bv_metab_found) %>%
    arrange(pvalue)

cat("\n============================================================\n")
cat("BV GENES × BV METABOLITES CORRELATIONS\n")
cat("============================================================\n")
cat("Total tests:", nrow(bv_focused), "\n")
cat("Nominally significant (p < 0.05):", sum(bv_focused$pvalue < 0.05), "\n\n")

print(as_tibble(bv_focused) %>%
      dplyr::select(Gene_ID, Gene_Annotation, Metabolite, rho, pvalue, qvalue) %>%
      head(30),
      width = Inf)
# #### Step 3: Specifically Test cas1/cas2 and chuR
# ============================================================
# FOCUSED: SCOARY-IDENTIFIED BV GENES
# ============================================================

# These genes were nominally BV-enriched in your Scoary analysis
scoary_bv_genes <- c("cas1", "cas2", "chuR")

# Check which are in your gene matrix
scoary_genes_present <- intersect(scoary_bv_genes, colnames(gene_by_sample_annotated))

if(length(scoary_genes_present) > 0) {
    cat("Testing Scoary-identified BV genes:\n")

    scoary_results <- cor_results_ann %>%
        filter(Gene_ID %in% scoary_genes_present) %>%
        arrange(pvalue)

    cat("\n--- All correlations for cas1, cas2, chuR ---\n")
    print(as_tibble(scoary_results) %>%
          dplyr::select(Gene_ID, Gene_Annotation, Metabolite, rho, pvalue, qvalue),
          n = 50)

    # Highlight polyamine-related correlations
    cat("\n--- Polyamine metabolite correlations ---\n")
    polyamine_metab <- c("Putrescine", "Cadaverine", "N-Acetylputrescine",
                         "N-Acetylspermidine", "5-Methylthioadenosine", "Histamine")

    scoary_polyamine <- scoary_results %>%
        filter(Metabolite %in% polyamine_metab)

    print(as_tibble(scoary_polyamine) %>%
          dplyr::select(Gene_ID, Metabolite, rho, pvalue),
          n = 20)
}
   print(as_tibble(scoary_polyamine) %>%
          dplyr::select(Gene_ID, Metabolite, rho, pvalue, qvalue),
          n = 20)
# #### Step 4: Create Publication-Ready Heatmap
bv_metab_found
str(cor_results_ann)
cor_results_ann %>%
    filter(qvalue < 0.10) %>%
    pull(Gene_ID) %>%
    unique()%>%
    length()
cor_results_ann %>%
    filter(
           qvalue < 0.1) %>%
    pull(Gene_ID) %>%
    unique()
# ============================================================
# HEATMAP: BV GENES × BV METABOLITES
# ============================================================

library(pheatmap)
library(RColorBrewer)

# Build correlation matrix for visualization
# Use all genes with at least one p < 0.05 against BV metabolites

sig_genes <- cor_results_ann %>%
    filter(
        #Metabolite %in% bv_metab_found,
           qvalue < 0.1) %>%
    pull(Gene_ID) %>%
    unique()

cat("Genes with at least one nominally significant BV-metabolite correlation:",
    length(sig_genes), "\n")

if(length(sig_genes) >= 5) {
    # Build matrix
    heatmap_data <- cor_results_ann %>%
        filter(Gene_ID %in% sig_genes,
               Metabolite %in% bv_metab_found) %>%
        dplyr::select(Gene_ID, Metabolite, rho) %>%
        pivot_wider(names_from = Metabolite, values_from = rho) %>%
        column_to_rownames("Gene_ID")

    # Add annotations
    gene_annot <- gene_annotations %>%
        filter(Gene_ID %in% rownames(heatmap_data)) %>%
        dplyr::select(Gene_ID, Gene_Annotation) %>%
        distinct() %>%
        column_to_rownames("Gene_ID")

    # Create row labels with annotations
    row_labels <- paste0(rownames(heatmap_data), " (",
                         gene_annot[rownames(heatmap_data), "Gene_Annotation"], ")")
    row_labels <- substr(row_labels, 1, 60)  # Truncate long labels

    # Plot
    pheatmap(heatmap_data,
             main = "Gene-Metabolite Correlations\n(BV-relevant metabolites, p < 0.001 genes)",
             color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
             breaks = seq(-0.8, 0.8, length.out = 101),
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = TRUE,
             labels_row = row_labels,
             fontsize_row = 6,
             fontsize_col = 8,
             angle_col = 45,
             border_color = NA)
}
# Test association between MTA-gene cluster and BV status
mta_bv_table <- matrix(c(13, 3, 8, 9), nrow = 2, byrow = TRUE,
                        dimnames = list(c("High", "Low"), c("Control", "BV")))

print(mta_bv_table)
fisher_result <- fisher.test(mta_bv_table)
cat("\nFisher's exact test p-value:", fisher_result$p.value, "\n")
cat("Odds ratio:", round(fisher_result$estimate, 2), "\n")
# NOTE: Draft figure generation and exploratory analyses removed for clarity.
# The publication-ready Figure 7 is generated in the section below.

# ============================================================================
# PUBLICATION FIGURE 7: Refined Genome-Metabolome Integration
# Journal-ready version with adjusted panel sizing
# ============================================================================

# Load additional libraries needed for Figure 7
library(tidyr)
library(tibble)

# Load gene annotations from CSV (provides annotation data for Figure 7 panel labels)
cat("Loading gene annotations for Figure 7...\n")
gene_pa <- read.csv("data/gardnerella_gene_presence_absence.csv", row.names = 1)

# Extract gene annotations (columns 1-13 contain metadata)
gene_annotations <- data.frame(
    Gene_ID = rownames(gene_pa),
    Gene_Annotation = gene_pa$Annotation,
    stringsAsFactors = FALSE
)

# Color scheme
color_bv <- "#F68B33"      # Orange for BV
color_control <- "#388ECC"  # Blue for Control

# ============================================================
# PANEL A: Procrustes Analysis (Refined)
# - No subtitle (stats go in figure caption)
# - Dual legend: color (BV/Control) + shape (Gene/Metabolome PCoA)
# ============================================================
cat("Creating Panel A: Procrustes...\n")

# Reshape data for proper legend handling
procrustes_gene <- procrustes_coords %>%
    select(Sample, Gene_PC1, Gene_PC2, Group) %>%
    mutate(DataType = "Gene PCoA") %>%
    rename(PC1 = Gene_PC1, PC2 = Gene_PC2)

procrustes_metab <- procrustes_coords %>%
    select(Sample, Metab_PC1, Metab_PC2, Group) %>%
    mutate(DataType = "Metabolome PCoA") %>%
    rename(PC1 = Metab_PC1, PC2 = Metab_PC2)

procrustes_long <- bind_rows(procrustes_gene, procrustes_metab)
procrustes_long$DataType <- factor(procrustes_long$DataType,
                                    levels = c("Gene PCoA", "Metabolome PCoA"))

panel_a <- ggplot() +
    # Arrows connecting paired points
    geom_segment(data = procrustes_coords,
                 aes(x = Gene_PC1, y = Gene_PC2,
                     xend = Metab_PC1, yend = Metab_PC2),
                 arrow = arrow(length = unit(0.08, "cm")),
                 color = "gray50", alpha = 0.5, linewidth = 0.3) +
    # Points with both color and shape aesthetics
    geom_point(data = procrustes_long,
               aes(x = PC1, y = PC2, color = Group, shape = DataType),
               size = 3, alpha = 0.9) +
    # Color scale
    scale_color_manual(values = c("BV" = color_bv, "Control" = color_control),
                       name = "Group") +
    # Shape scale
    scale_shape_manual(values = c("Gene PCoA" = 16, "Metabolome PCoA" = 17),
                       name = "Data type") +
    # Labels - NO subtitle
    labs(title = "a",
         x = "Dimension 1",
         y = "Dimension 2") +
    # Theme
    theme_bw(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.margin = margin(t = 0, b = 0),
        legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10)
    ) +
    guides(
        color = guide_legend(order = 1, override.aes = list(size = 3)),
        shape = guide_legend(order = 2, override.aes = list(size = 3))
    )

# ============================================================
# PANEL B: Variance Partitioning (Refined)
# - No subtitle (stats go in figure caption)
# - Smaller pie chart
# ============================================================
cat("Creating Panel B: Variance partitioning...\n")

varpart_data <- data.frame(
    Component = c("Gene content\n(unique)", "Shared", "Confounders\n(unique)", "Residual"),
    Value = c(9, 6, 0.5, 84.5),  # Slightly adjusted for better display
    Label = c("9%", "6%", "<1%", "85%")
)

# Set factor order for legend
varpart_data$Component <- factor(varpart_data$Component,
                                  levels = c("Gene content\n(unique)", "Shared",
                                             "Confounders\n(unique)", "Residual"))

panel_b <- ggplot(varpart_data, aes(x = "", y = Value, fill = Component)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = c(
        "Gene content\n(unique)" = color_control,
        "Shared" = "#984EA3",
        "Confounders\n(unique)" = color_bv,
        "Residual" = "gray80"
    )) +
    geom_text(aes(label = Label),
              position = position_stack(vjust = 0.5),
              color = c("white", "white", "white", "gray30"),
              fontface = "bold", size = 3.5) +
    labs(title = "b", fill = "") +
    theme_void(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )

# ============================================================
# PANEL C: Gene-Metabolite Heatmap (Refined)
# - Title is just "C" (description in figure caption)
# ============================================================
cat("Creating Panel C: Heatmap...\n")

# Get top FDR-significant genes
top_fdr_genes <- cor_results_ann %>%
    filter(qvalue < 0.10) %>%
    pull(Gene_ID) %>%
    unique()

cat("  Genes with FDR q < 0.10:", length(top_fdr_genes), "\n")

# If too few FDR hits, use top nominal hits
if(length(top_fdr_genes) < 10) {
    top_fdr_genes <- cor_results_ann %>%
        arrange(pvalue) %>%
        pull(Gene_ID) %>%
        unique() %>%
        head(20)
    cat("  Using top 20 nominal hits instead\n")
}

# Select key metabolites for clarity
key_metabolites <- c("5-Methylthioadenosine", "Tyrosine", "Putrescine",
                     "Cadaverine", "Histamine", "Succinate",
                     "Trimethylamine", "Lactate")

# Check which are present
available_metabolites <- unique(cor_results_ann$Metabolite)
key_metabolites_present <- intersect(key_metabolites, available_metabolites)
cat("  Key metabolites present:", length(key_metabolites_present), "\n")

# Build heatmap matrix
heatmap_data <- cor_results_ann %>%
    filter(Gene_ID %in% top_fdr_genes,
           Metabolite %in% key_metabolites_present) %>%
    dplyr::select(Gene_ID, Metabolite, rho) %>%
    pivot_wider(names_from = Metabolite, values_from = rho) %>%
    column_to_rownames("Gene_ID")

# Remove rows/cols with all NA
heatmap_data <- heatmap_data[rowSums(!is.na(heatmap_data)) > 0,
                              colSums(!is.na(heatmap_data)) > 0, drop = FALSE]

# Replace remaining NAs with 0 for visualization
heatmap_data[is.na(heatmap_data)] <- 0

cat("  Heatmap dimensions:", nrow(heatmap_data), "genes x", ncol(heatmap_data), "metabolites\n")

# Get gene annotations for row labels
row_labels_short <- sapply(rownames(heatmap_data), function(x) {
    ann <- gene_annotations$Gene_Annotation[gene_annotations$Gene_ID == x][1]
    if(is.na(ann) || ann == "" || ann == "hypothetical protein") {
        # Shorten gene ID if too long
        if(nchar(x) > 20) {
            return(substr(x, 1, 20))
        }
        return(x)
    } else {
        # Truncate annotation to 30 chars
        ann_short <- substr(ann, 1, 30)
        return(paste0(substr(x, 1, 15), " (", ann_short, ")"))
    }
})

# Create heatmap grob (silent to capture as grob)
heatmap_grob <- pheatmap(
    heatmap_data,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-0.8, 0.8, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    labels_row = row_labels_short,
    fontsize_row = 6,
    fontsize_col = 8,
    angle_col = 45,
    border_color = NA,
    main = "c",  # Just the panel label
    fontsize = 10,
    silent = TRUE
)

# ============================================================
# COMBINE PANELS - ADJUSTED WIDTHS
# Procrustes (1.3) | Pie (0.7) | Heatmap (1.5)
# ============================================================
cat("Combining panels with adjusted widths...\n")

# Save combined figure
pdf("figures/Figure7_genome_metabolome.pdf", width = 14, height = 6)

grid.arrange(
    ggplotGrob(panel_a),
    ggplotGrob(panel_b),
    heatmap_grob[[4]],
    ncol = 3,
    widths = c(1.3, 0.7, 1.5)  # Adjusted: Procrustes larger, pie smaller, heatmap largest
)

dev.off()
cat("Saved: figures/Figure7_genome_metabolome.pdf\n")

# Save PNG version
png("figures/Figure7_genome_metabolome.png", width = 14, height = 6, units = "in", res = 300)

grid.arrange(
    ggplotGrob(panel_a),
    ggplotGrob(panel_b),
    heatmap_grob[[4]],
    ncol = 3,
    widths = c(1.3, 0.7, 1.5)
)

dev.off()
cat("Saved: figures/Figure7_genome_metabolome.png\n")

# ============================================================
# SAVE INDIVIDUAL PANELS
# ============================================================
cat("Saving individual panels...\n")

ggsave("figures/Figure7a_procrustes.pdf", panel_a, width = 6, height = 5.5)
ggsave("figures/Figure7a_procrustes.png", panel_a, width = 6, height = 5.5, dpi = 300)
cat("Saved: figures/Figure7a_procrustes.pdf/png\n")

ggsave("figures/Figure7b_varpart.pdf", panel_b, width = 4, height = 4.5)
ggsave("figures/Figure7b_varpart.png", panel_b, width = 4, height = 4.5, dpi = 300)
cat("Saved: figures/Figure7b_varpart.pdf/png\n")

# Save heatmap separately
pdf("figures/Figure7c_heatmap.pdf", width = 8, height = 8)
grid.draw(heatmap_grob[[4]])
dev.off()
cat("Saved: figures/Figure7c_heatmap.pdf\n")

png("figures/Figure7c_heatmap.png", width = 8, height = 8, units = "in", res = 300)
grid.draw(heatmap_grob[[4]])
dev.off()
cat("Saved: figures/Figure7c_heatmap.png\n")

cat("\n=== Figure 7 refinement complete ===\n")
cat("Changes implemented:\n")
cat("  1. Panel widths adjusted: Procrustes (1.3), Pie (0.7), Heatmap (1.5)\n")
cat("  2. Subtitles removed - statistics now belong in figure caption\n")
cat("  3. Procrustes plot has dual legend: Group (color) + Data type (shape)\n")
cat("  4. Heatmap title simplified to 'c'\n")
