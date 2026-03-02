# ===================================================================
# MitCH Study: MOFA2 Multi-Omics Factor Analysis
# ===================================================================
#
# Purpose:
#   Unsupervised multi-omics integration using MOFA2 with 4 data views
#   (metabolomics_pos, metabolomics_neg, pathways, taxonomy) for
#   comprehensive vaginal microbiome analysis.
#
# Inputs:
#   - objects/multiomics_views.rda      (X: list of 4 omics matrices)
#   - metadata/sample_metadata.rda      (meta_data: clinical metadata)
#   - objects/group_labels.rda          (Y_aligned: BV/Control labels)
#   - objects/diablo_selected_features.rda (all_selected_features from DIABLO)
#
# Outputs:
#   - figures/MOFA2_Factor1_loadings.pdf/png       (main text figure)
#   - figures/MOFA2_loadings_supplementary.pdf/png  (supplementary figure)
#   - figures/MOFA2_variance_explained.png
#   - figures/MOFA2_pH_analysis.png
#   - figures/MOFA2_metabolomics_pos.png
#   - figures/MOFA2_metabolomics_neg.png
#   - figures/MOFA2_pathways.png
#   - figures/MOFA2_taxonomy.png
#   - supplementary_tables/Table_S9.xlsx  (variance explained)
#   - supplementary_tables/Table_S10.xlsx (complete feature loadings)
#
# Key Features:
#   1. Your exact data preprocessing pipeline
#   2. Enhanced 4-view MOFA2 analysis
#   3. Metabolomics pathway interpretation
#   4. Clinical-microbiome integration
#   5. Advanced biomarker discovery
# ===================================================================

# Load required libraries
library(MOFA2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(corrplot)
library(RColorBrewer)
library(cowplot)
library(purrr)
library(ggtext)    # For italicized species names in plots
library(openxlsx)  # For Excel export of supplementary tables


# ============================================================================
# DATA PREPARATION WITH SEPARATED METABOLOMICS MODES (YOUR PIPELINE)
# ============================================================================

# Load your data with 4 views: metabolomics_pos, metabolomics_neg, pathways, taxonomy
load("objects/multiomics_views.rda")  # Assumes X contains the 4 views
load("metadata/sample_metadata.rda")
load("objects/group_labels.rda")

# Verify data structure
cat("Data structure overview:\n")
str(X)
glimpse(meta_data)




# Ensure sample order matches between omics data and metadata

X$metabolomics_pos <- X$metabolomics_pos[meta_data$Sample,]
X$metabolomics_neg <- X$metabolomics_neg[meta_data$Sample,]
X$pathways <- X$pathways[meta_data$Sample,]
X$taxonomy <- X$taxonomy[meta_data$Sample,]
rowSums(X$pathways)
rowSums(X$taxonomy)
rowSums(X$metabolomics_neg)

sample_order <- rownames(X$metabolomics_pos)  # Use pos mode as reference
meta_data_aligned <- meta_data[match(sample_order, meta_data$Sample), ]

# Verify alignment
if(!all(meta_data_aligned$Sample == sample_order)) {
  stop("Sample order mismatch between omics data and metadata!")
}

print("Sample alignment verified. Metadata summary:")
print(table(meta_data_aligned$Group))
print(summary(meta_data_aligned[, c("Age", "Vaginal_pH")]))

# Create comprehensive sample metadata (ENHANCED)
sample_metadata <- meta_data_aligned


rownames(sample_metadata) <- sample_metadata$Sample
glimpse(sample_metadata)

# ============================================================================
# ENHANCED PREPROCESSING FOR 4-VIEW DATA
# ============================================================================

# Transform compositional data (taxonomy and pathways)
# Using CLR transformation for compositional data
clr_transform <- function(x) {
  # Add pseudocount for zeros
  x_pseudo <- x + min(x[x > 0]) * 0.01
  # Compute geometric mean per sample
  gm <- exp(rowMeans(log(x_pseudo)))
  # CLR transformation
  log(x_pseudo / gm)
}

# Prepare data blocks
#X_prepared <- list(
#  metabolomics_pos =  t(scale(t(X$metabolomics_pos), center = TRUE, scale = TRUE)),
#  metabolomics_neg =  t(scale(t(X$metabolomics_neg), center = TRUE, scale = TRUE)),
#  pathways = t(scale(t(X$pathways), center = TRUE, scale = TRUE)),
#  taxonomy = t(scale(t(X$taxonomy), center = TRUE, scale = TRUE))
#)
# Prepare data blocks
X_prepared <- list(
  metabolomics_pos =  X$metabolomics_pos,
  metabolomics_neg =  X$metabolomics_neg,
  pathways = X$pathways,
  taxonomy = X$taxonomy
)


# Check for NAs and infinite values
check_data_quality <- function(data_list) {
  lapply(names(data_list), function(name) {
    x <- data_list[[name]]
    tibble::tibble(
      Block = name,
      N_samples = nrow(x),
      N_features = ncol(x),
      N_NA = sum(is.na(x)),
      N_Inf = sum(is.infinite(x)),
      Min = min(x, na.rm = TRUE),
      Max = max(x, na.rm = TRUE)
    )
  }) %>% dplyr::bind_rows()
}

data_quality <- check_data_quality(X_prepared)
print("Data Quality Check:")
print(data_quality)

# Replace any remaining NAs/Infs with 0
X_prepared <- lapply(X_prepared, function(x) {
  x[is.na(x) | is.infinite(x)] <- 0
  x
})

X_processed <- X_prepared

# ============================================================================
# MOFA2 MODEL SETUP WITH 4 VIEWS (YOUR PIPELINE + ENHANCEMENTS)
# ============================================================================

# Ensure all views have the same samples in the same order
common_samples <- rownames(X_processed$metabolomics_pos)
for(view in names(X_processed)) {
  common_samples <- intersect(common_samples, rownames(X_processed[[view]]))
}
cat("Common samples across all views:", length(common_samples), "\n")

# Subset all views to common samples and ensure same order
X_aligned <- list()
for(view in names(X_processed)) {
  X_aligned[[view]] <- X_processed[[view]][common_samples, ]
  cat(view, "dimensions after alignment:", paste(dim(X_aligned[[view]]), collapse = " x "), "\n")
}

# Transpose matrices for MOFA2 (features as rows, samples as columns)
X_mofa <- list()
for(view in names(X_aligned)) {
  X_mofa[[view]] <- t(X_aligned[[view]])
  cat(view, "MOFA format:", paste(dim(X_mofa[[view]]), collapse = " x "),
      "(features x samples)\n")
}

# Verify sample order consistency
sample_orders_match <- TRUE
for(i in 2:length(X_mofa)) {
  if(!identical(colnames(X_mofa[[1]]), colnames(X_mofa[[i]]))) {
    sample_orders_match <- FALSE
    break
  }
}

if(!sample_orders_match) {
  stop("Sample orders don't match across views!")
} else {
  cat("Sample order verified across all views.\n")
}

# Create MOFA object
MOFAobject <- create_mofa(X_mofa)

# Align metadata to common samples
sample_metadata_aligned <- sample_metadata[common_samples, ]
mofa_samples <- colnames(X_mofa[[1]])
row.names(sample_metadata_aligned) <- sample_metadata_aligned$Sample
sample_metadata_mofa <- sample_metadata_aligned[mofa_samples, ]

# Add required 'sample' column for MOFA2
sample_metadata_mofa$sample <- mofa_samples

# Verify perfect alignment
if(!identical(sample_metadata_mofa$sample, colnames(X_mofa[[1]]))) {
  stop("Sample alignment failed!")
} else {
  cat("Sample alignment verified: metadata and MOFA samples match perfectly.\n")
}



# Add comprehensive metadata
samples_metadata(MOFAobject) <- sample_metadata_mofa
print("Data overview:")
print(MOFAobject)

# ============================================================================
# ENHANCED MODEL OPTIONS FOR 4-VIEW ANALYSIS
# ============================================================================

# Define data options
data_opts <- get_default_data_options(MOFAobject)

# Scale views independently (important for different data types)
data_opts$scale_views <- FALSE

# Don't center groups (single group study)
#data_opts$center_groups <- FALSE

# Define model options with 4 views
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10  # Increase for more complex patterns with 4 views
model_opts$likelihoods <- c(
  "metabolomics_pos" = "gaussian",
  "metabolomics_neg" = "gaussian",
  "pathways" = "gaussian",
  "taxonomy" = "gaussian"
)
# OPTIMAL sparsity settings for biological interpretation
# Enable sparse factors for cleaner interpretation
#model_opts$spikeslab_factors <- FALSE  # Changed from FALSE

# Enable sparse weights for feature selection
model_opts$spikeslab_weights <- TRUE

# ARD (Automatic Relevance Determination) hyperparameters
# These control sparsity strength
#model_opts$ard_factors <- TRUE
#model_opts$ard_weights <- TRUE

# Enhanced training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 123
train_opts$freqELBO <- 1
train_opts$maxiter <- 3000


# Early stopping parameters
train_opts$startELBO <- 100  # Start checking after 100 iterations
#train_opts$tolerance <- 0.01  # Stop if improvement < 1%

# OPTIMAL factor dropping threshold
# Higher threshold = more aggressive factor removal
train_opts$drop_factor_threshold <- 0.02  # Increased from 0.015

# GPU settings (if available)
train_opts$gpu_mode <- FALSE  # Set TRUE if CUDA available

# Stochastic settings for large datasets (not needed here)
train_opts$stochastic <- FALSE

# Train the model
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)



print("Training MOFA2 model with 4 omics views...")
MOFAobject <- run_mofa(MOFAobject, use_basilisk = TRUE)



library(caret)
library(randomForest)
# ============================================================================
# ENHANCED ANALYSIS FUNCTIONS FOR 4-VIEW DATA
# ============================================================================

# Enhanced metabolite interpretation function
interpret_metabolite_biological_role <- function(metabolite_name, ion_mode) {

  metabolite_lower <- tolower(metabolite_name)

  # Amino acids and derivatives
  if(grepl("acetylcarnitine|carnitine", metabolite_lower)) {
    if(grepl("c2:|acetyl", metabolite_lower)) {
      return("Short-chain acylcarnitine - acetyl-CoA metabolism, energy production")
    } else {
      return("Acylcarnitine - fatty acid \u03b2-oxidation, mitochondrial function")
    }
  }
  else if(grepl("alanine|leucine|isoleucine|valine|tryptophan|tyrosine|phenylalanine|methionine", metabolite_lower)) {
    if(grepl("acetyl", metabolite_lower)) {
      return("N-acetyl amino acid - amino acid metabolism, detoxification")
    } else {
      return("Amino acid - protein metabolism, neurotransmitter precursor")
    }
  }
  else if(grepl("gaba|glutam", metabolite_lower)) {
    return("Neurotransmitter/neuromodulator - GABA pathway, stress response")
  }

  # Nucleotides and nucleosides
  else if(grepl("cmp|ump|gmp|amp|cytidine|uridine|adenosine|guanosine", metabolite_lower)) {
    return("Nucleotide/nucleoside - DNA/RNA synthesis, energy metabolism")
  }

  # Organic acids (important for vaginal environment)
  else if(grepl("lactate|succinate|fumarate|malate|citrate", metabolite_lower)) {
    if(grepl("lactate", metabolite_lower)) {
      return("Lactate - lactobacilli fermentation, vaginal acidification")
    } else {
      return("TCA cycle intermediate - central energy metabolism")
    }
  }
  else if(grepl("aconitate", metabolite_lower)) {
    return("TCA cycle intermediate - energy metabolism")
  }

  # Host waste products and detoxification
  else if(grepl("allantoin", metabolite_lower)) {
    return("Purine catabolite - oxidative stress marker, antioxidant")
  }
  else if(grepl("cysteine", metabolite_lower)) {
    return("Sulfur amino acid - antioxidant defense, detoxification")
  }

  # Ion mode specific interpretations
  mode_note <- ifelse(ion_mode == "metabolomics_pos",
                      " [ESI+: basic/polar compounds]",
                      " [ESI-: acidic/anionic compounds]")

  return(paste0("Metabolite", mode_note))
}

# Enhanced factor analysis for 4 views
analyze_enhanced_factor_loadings <- function(MOFAobject, factor_num = 1, top_n = 20) {

  cat(sprintf("\n=== ENHANCED FACTOR %d ANALYSIS - 4 OMICS VIEWS ===\n", factor_num))

  weights <- get_weights(MOFAobject)
  all_loadings <- data.frame()

  for(view_name in names(weights)) {
    view_weights <- weights[[view_name]][, factor_num]

    view_df <- data.frame(
      Feature = names(view_weights),
      Loading = view_weights,
      Abs_Loading = abs(view_weights),
      View = view_name,
      stringsAsFactors = FALSE
    )

    all_loadings <- rbind(all_loadings, view_df)
  }

  # Sort by absolute loading strength
  all_loadings <- all_loadings[order(-all_loadings$Abs_Loading), ]
  top_features <- head(all_loadings, top_n)

  cat(sprintf("Top %d features driving Factor %d:\n", top_n, factor_num))
  print(top_features[, c("Feature", "Loading", "View")])

  # Enhanced interpretation by view
  cat(sprintf("\n=== BIOLOGICAL INTERPRETATION BY OMICS VIEW ===\n"))

  for(view in unique(all_loadings$View)) {
    view_features <- all_loadings[all_loadings$View == view, ]
    cat(sprintf("\n%s:\n", toupper(view)))

    top_positive <- head(view_features[view_features$Loading > 0, ], 5)
    top_negative <- head(view_features[view_features$Loading < 0, ], 5)

    if(nrow(top_positive) > 0) {
      cat("  UP POSITIVE LOADINGS (higher in high factor score samples):\n")
      for(i in 1:nrow(top_positive)) {
        feature_name <- top_positive$Feature[i]
        loading_val <- top_positive$Loading[i]

        if(view %in% c("metabolomics_pos", "metabolomics_neg")) {
          interpretation <- interpret_metabolite_biological_role(feature_name, view)
          cat(sprintf("    * %s (%.3f) - %s\n", feature_name, loading_val, interpretation))
        } else {
          cat(sprintf("    * %s (%.3f)\n", feature_name, loading_val))
        }
      }
    }

    if(nrow(top_negative) > 0) {
      cat("  DOWN NEGATIVE LOADINGS (lower in high factor score samples):\n")
      for(i in 1:nrow(top_negative)) {
        feature_name <- top_negative$Feature[i]
        loading_val <- top_negative$Loading[i]

        if(view %in% c("metabolomics_pos", "metabolomics_neg")) {
          interpretation <- interpret_metabolite_biological_role(feature_name, view)
          cat(sprintf("    * %s (%.3f) - %s\n", feature_name, loading_val, interpretation))
        } else {
          cat(sprintf("    * %s (%.3f)\n", feature_name, loading_val))
        }
      }
    }
  }

  return(all_loadings)
}

# Variance explained analysis for 4 views
analyze_4view_variance_explained <- function(MOFAobject) {

  cat("\n=== 4-VIEW VARIANCE EXPLAINED ANALYSIS ===\n")

  variance_explained <- get_variance_explained(MOFAobject)$r2_per_factor$group1

  cat("Variance explained by each factor across all views:\n")
  print(round(variance_explained, 2))

  # Calculate total variance explained per view
  total_variance <- apply(variance_explained, 1, sum)
  cat("\nTotal variance explained per view:\n")
  for(i in 1:length(total_variance)) {
    cat(sprintf("* %s: %.2f%%\n", names(total_variance)[i], total_variance[i]))
  }

  # Create enhanced visualization
  variance_df <- as.data.frame(variance_explained)
  variance_df$Factor <- rownames(variance_df)
  variance_long <- tidyr::gather(variance_df, View, Variance, -Factor)

  variance_df
  # Enhanced color scheme for 4 views
  view_colors <- c(
    "metabolomics_pos" = "#E91E63",  # Pink
    "metabolomics_neg" = "#3F51B5",  # Indigo
    "pathways" = "#4CAF50",          # Green
    "taxonomy" = "#FF9800"           # Orange
  )

  p_variance <- ggplot(variance_long, aes(x = Factor, y = Variance, fill = View)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = view_colors) +
    labs(title = "Variance Explained by MOFA Factors Across 4 Omics Views",
         subtitle = "Metabolomics (Pos/Neg), Functional Pathways, and Taxonomy",
         x = "MOFA Factors",
         y = "Variance Explained (%)",
         fill = "View") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0),
          legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 1))

  print(p_variance)

  ggsave("figures/MOFA2_variance_explained.png", p_variance, width = 8, height = 6)
  return(list(variance_table = variance_explained, variance_plot = p_variance))
}

# Enhanced microbiome state classification using 4 views
classify_enhanced_microbiome_states <- function(MOFAobject, sample_metadata) {

  cat("\n=== ENHANCED MICROBIOME STATE CLASSIFICATION ===\n")
  cat("Using 4-view MOFA factors for biological state definition\n")

  factors <- get_factors(MOFAobject)[[1]]

  # Initialize enhanced microbiome state
  sample_metadata$Enhanced_Microbiome_State <- "Unknown"

  # Enhanced classification using multiple factors
  for(i in 1:nrow(factors)) {
    factor1 <- factors[i, "Factor1"]
    factor2 <- factors[i, "Factor2"]
    factor3 <- if("Factor3" %in% colnames(factors)) factors[i, "Factor3"] else 0

    # More sophisticated state classification
    if(factor1 < -0.8) {
      sample_metadata$Enhanced_Microbiome_State[i] <- "Optimal_Lactobacilli"
    } else if(factor1 < -0.3) {
      sample_metadata$Enhanced_Microbiome_State[i] <- "Healthy_Lactobacilli"
    } else if(factor2 > 0.5 && factor1 < 0.3) {
      sample_metadata$Enhanced_Microbiome_State[i] <- "Transitional_Iners"
    } else if(factor1 > 0.5) {
      if(factor3 > 0.3) {
        sample_metadata$Enhanced_Microbiome_State[i] <- "Severe_BV"
      } else {
        sample_metadata$Enhanced_Microbiome_State[i] <- "Moderate_BV"
      }
    } else {
      sample_metadata$Enhanced_Microbiome_State[i] <- "Mixed_State"
    }
  }

  # Convert to ordered factor
  sample_metadata$Enhanced_Microbiome_State <- factor(
    sample_metadata$Enhanced_Microbiome_State,
    levels = c("Optimal_Lactobacilli", "Healthy_Lactobacilli",
               "Mixed_State", "Transitional_Iners", "Moderate_BV", "Severe_BV")
  )

  cat("Enhanced microbiome state distribution:\n")
  print(table(sample_metadata$Enhanced_Microbiome_State))

  return(sample_metadata)
}

# Enhanced classification analysis with 4 views
run_enhanced_4view_classification <- function(MOFAobject, sample_metadata) {

  cat("\n=== ENHANCED 4-VIEW CLASSIFICATION ANALYSIS ===\n")

  # Extract MOFA factors
  mofa_factors <- get_factors(MOFAobject)[[1]]

  # Prepare enhanced prediction data
  prediction_data <- data.frame(
    mofa_factors,
    Age = as.numeric(sample_metadata$Age),
    Vaginal_pH = as.numeric(sample_metadata$Vaginal_pH),
    BV_History = as.numeric(sample_metadata$BV_History),
    Diagnosis = sample_metadata$Group
  )

  # Remove incomplete cases
  complete_cases <- complete.cases(prediction_data)
  prediction_data_clean <- prediction_data[complete_cases, ]
  prediction_data_clean$Diagnosis <- droplevels(prediction_data_clean$Diagnosis)

  cat(sprintf("Classification dataset: %d samples, %d classes\n",
              nrow(prediction_data_clean), nlevels(prediction_data_clean$Diagnosis)))
  print(table(prediction_data_clean$Diagnosis))

  results <- list()

  # Strategy 1: 4-View MOFA Factors Only
  cat("\n--- 4-View MOFA Factors Classification ---\n")
  factors_data <- prediction_data_clean[, c(grep("Factor", names(prediction_data_clean)),
                                            which(names(prediction_data_clean) == "Diagnosis"))]

  set.seed(123)
  rf_4view <- randomForest(Diagnosis ~ .,
                           data = factors_data,
                           ntree = 1000,
                           importance = TRUE,
                           proximity = TRUE)

  results$rf_4view <- rf_4view
  cat(sprintf("4-View MOFA OOB Error: %.3f (Accuracy: %.3f)\n",
              rf_4view$err.rate[nrow(rf_4view$err.rate), 1],
              1 - rf_4view$err.rate[nrow(rf_4view$err.rate), 1]))

  # Strategy 2: View-specific analysis
  weights <- get_weights(MOFAobject)

  for(view_name in names(weights)) {
    cat(sprintf("\n--- %s-Specific Classification ---\n", toupper(view_name)))

    # Get top features from this view
    view_weights <- weights[[view_name]]
    max_loadings <- apply(abs(view_weights), 1, max)
    top_features_idx <- order(max_loadings, decreasing = TRUE)[1:min(20, length(max_loadings))]

    # Create view-specific factor scores
    view_factors <- data.frame(mofa_factors)
    names(view_factors) <- paste0(view_name, "_", names(view_factors))

    view_data <- cbind(view_factors, Diagnosis = prediction_data_clean$Diagnosis)

    set.seed(123)
    rf_view <- randomForest(Diagnosis ~ .,
                            data = view_data,
                            ntree = 1000,
                            importance = TRUE)

    results[[paste0("rf_", view_name)]] <- rf_view
    cat(sprintf("%s OOB Error: %.3f (Accuracy: %.3f)\n",
                view_name,
                rf_view$err.rate[nrow(rf_view$err.rate), 1],
                1 - rf_view$err.rate[nrow(rf_view$err.rate), 1]))
  }

  # Strategy 3: Combined with clinical
  cat("\n--- Combined 4-View + Clinical ---\n")
  set.seed(123)
  rf_combined <- randomForest(Diagnosis ~ .,
                              data = prediction_data_clean,
                              ntree = 1000,
                              importance = TRUE)

  results$rf_combined <- rf_combined
  cat(sprintf("Combined OOB Error: %.3f (Accuracy: %.3f)\n",
              rf_combined$err.rate[nrow(rf_combined$err.rate), 1],
              1 - rf_combined$err.rate[nrow(rf_combined$err.rate), 1]))

  return(results)
}

# Enhanced visualization for 4 views
create_enhanced_4view_heatmaps <- function(MOFAobject, top_n = 20) {

  cat("\n=== CREATING ENHANCED 4-VIEW HEATMAPS ===\n")

  weights <- get_weights(MOFAobject)
  heatmap_list <- list()

  # Define color schemes for each view
  view_colors <- list(
    metabolomics_pos = colorRampPalette(c("#1565C0", "white", "#E91E63"))(100),
    metabolomics_neg = colorRampPalette(c("#2E7D32", "white", "#E65100"))(100),
    pathways = colorRampPalette(c("#4A148C", "white", "#FF6F00"))(100),
    taxonomy = colorRampPalette(c("#B71C1C", "white", "#1B5E20"))(100)
  )

  for(view_name in names(weights)) {

    cat(sprintf("Creating heatmap for %s...\n", view_name))

    view_matrix <- weights[[view_name]]

    if(nrow(view_matrix) == 0) {
      cat(sprintf("No features in %s\n", view_name))
      next
    }

    # Get top features by variance across factors
    max_loadings <- apply(abs(view_matrix), 1, max)
    top_idx <- order(max_loadings, decreasing = TRUE)[1:min(top_n, length(max_loadings))]
    heatmap_matrix <- view_matrix[top_idx, , drop = FALSE]

    # Enhanced annotation for biological context
    feature_names <- rownames(heatmap_matrix)

    if(view_name %in% c("metabolomics_pos", "metabolomics_neg")) {
      # Metabolite categories
      feature_categories <- sapply(feature_names, function(name) {
        name_lower <- tolower(name)
        if(grepl("carnitine|acetylcarnitine", name_lower)) {
          return("Fatty_Acid_Metabolism")
        } else if(grepl("amino|alanine|leucine|methionine|phenylalanine", name_lower)) {
          return("Amino_Acid_Metabolism")
        } else if(grepl("nucleotide|cmp|ump|gmp|amp|cytidine|uridine", name_lower)) {
          return("Nucleotide_Metabolism")
        } else if(grepl("lactate|citrate|succinate|aconitate", name_lower)) {
          return("Central_Carbon_Metabolism")
        } else {
          return("Other_Metabolites")
        }
      })
    } else if(view_name == "pathways") {
      # Pathway categories
      feature_categories <- sapply(feature_names, function(name) {
        name_lower <- tolower(name)
        if(grepl("ferment|glycol", name_lower)) {
          return("Fermentation_Pathways")
        } else if(grepl("amino|arg|met", name_lower)) {
          return("Amino_Acid_Pathways")
        } else if(grepl("fold|vitamin", name_lower)) {
          return("Vitamin_Cofactor_Pathways")
        } else {
          return("Other_Pathways")
        }
      })
    } else {
      # Taxonomy categories
      feature_categories <- sapply(feature_names, function(name) {
        if(grepl("Lactobacillus", name)) {
          if(grepl("iners", name)) {
            return("Transition_Lactobacilli")
          } else {
            return("Protective_Lactobacilli")
          }
        } else if(grepl("Gardnerella|Fannyhessea|Atopobium", name)) {
          return("Classic_BV_Pathogens")
        } else if(grepl("Prevotella|Sneathia|Dialister|Megasphaera", name)) {
          return("BV_Associated_Anaerobes")
        } else {
          return("Other_Taxa")
        }
      })
    }

    # Create annotation
    row_annotation <- data.frame(Category = feature_categories,
                                 row.names = feature_names)

    # Create heatmap
    pheat <- pheatmap(heatmap_matrix,
             annotation_row = row_annotation,
             scale = "none",
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = view_colors[[view_name]],
             main = paste0("Top ", nrow(heatmap_matrix), " Features - ",
                           toupper(view_name)),
             fontsize_row = 6,
             fontsize_col = 10)


    png(file = paste0("figures/MOFA2_", view_name, ".png"), width = 8, height = 7, units = "in", res = 300)
    print(pheat)
    dev.off()

    heatmap_list[[view_name]] <- heatmap_matrix
  }

  return(heatmap_list)
}

# ============================================================================
# MOFA LOADINGS DIVERGING LOLLIPOP FIGURE
# ============================================================================
# Creates publication-ready diverging lollipop plot showing top MOFA factor
# loadings across all 4 omics views, with DIABLO overlap annotation.
#
# Output:
#   - Main text figure: Factor 1 only (1 row x 4 cols)
#   - Supplementary figure: Factor 1 + Factor 2 (2 rows x 4 cols)
# ============================================================================

create_mofa_loadings_figure <- function(
    MOFAobject,
    factors = c(1),
    top_n = 5,
    diablo_features = NULL,
    figure_type = c("main", "supplementary"),
    output_dir = "figures"
) {

  figure_type <- match.arg(figure_type)

  cat("\n=== CREATING MOFA LOADINGS DIVERGING LOLLIPOP FIGURE ===\n")

  cat(sprintf("Figure type: %s | Factors: %s | Top N: %d\n",
              figure_type, paste(factors, collapse = ", "), top_n))

 # --------------------------------------------------------------------------
  # 1. Extract weights and variance explained
  # --------------------------------------------------------------------------
  weights <- MOFA2::get_weights(MOFAobject)
  var_explained <- MOFA2::get_variance_explained(MOFAobject)$r2_per_factor$group1

  # View display names (prettier labels for facets)
  view_labels <- c(
    "taxonomy" = "Taxonomy",
    "metabolomics_pos" = "Metabolites (ESI+)",
    "metabolomics_neg" = "Metabolites (ESI\u2212)",
    "pathways" = "Pathways"
  )

  # --------------------------------------------------------------------------
  # 2. Build long-format data frame with top features per factor x view
  # --------------------------------------------------------------------------
  all_loadings <- purrr::map_dfr(factors, function(f) {
    factor_name <- paste0("Factor", f)

    purrr::map_dfr(names(weights), function(view_name) {
      view_weights <- weights[[view_name]][, f]

      df <- dplyr::tibble(
        Feature = names(view_weights),
        Loading = as.numeric(view_weights),
        View = view_name,
        Factor = factor_name
      ) %>%
        dplyr::mutate(
          Abs_Loading = abs(Loading),
          Direction = dplyr::if_else(Loading > 0, "BV-enriched", "Control-enriched")
        )

      # Select top N positive + top N negative
      top_pos <- df %>%
        dplyr::filter(Loading > 0) %>%
        dplyr::arrange(dplyr::desc(Abs_Loading)) %>%
        dplyr::slice_head(n = top_n)

      top_neg <- df %>%
        dplyr::filter(Loading < 0) %>%
        dplyr::arrange(dplyr::desc(Abs_Loading)) %>%
        dplyr::slice_head(n = top_n)

      dplyr::bind_rows(top_pos, top_neg)
    })
  })

  # --------------------------------------------------------------------------
  # 3. Add DIABLO overlap indicator (no asterisks in labels)
  # --------------------------------------------------------------------------
  if (!is.null(diablo_features) && length(diablo_features) > 0) {
    all_loadings <- all_loadings %>%
      dplyr::mutate(
        DIABLO_overlap = Feature %in% diablo_features,
        Feature_label = Feature  # No asterisks
      )
  } else {
    all_loadings <- all_loadings %>%
      dplyr::mutate(
        DIABLO_overlap = FALSE,
        Feature_label = Feature
      )
  }

  # --------------------------------------------------------------------------
  # 4. Clean feature names for display
  # --------------------------------------------------------------------------
  # Helper function to format feature names
  format_feature_name <- function(name, view) {
    # For taxonomy: italicize species names using ggtext markdown
    if (view == "taxonomy") {
      # Extract genus and species, italicize
      if (grepl("^[A-Z][a-z]+\\s", name)) {
        # Standard binomial: Lactobacillus crispatus -> *L. crispatus*
        parts <- strsplit(name, " ")[[1]]
        genus_abbrev <- paste0(substr(parts[1], 1, 1), ".")
        species <- paste(parts[-1], collapse = " ")
        return(paste0("*", genus_abbrev, " ", species, "*"))
      } else if (grepl("^[A-Z][a-z]+_", name)) {
        # Underscore format: Lactobacillus_crispatus
        parts <- strsplit(name, "_")[[1]]
        genus_abbrev <- paste0(substr(parts[1], 1, 1), ".")
        species <- paste(parts[-1], collapse = " ")
        return(paste0("*", genus_abbrev, " ", species, "*"))
      }
    }

    # For pathways: abbreviate long names
    if (view == "pathways") {
      if (nchar(name) > 35) {
        # Truncate and add ellipsis
        return(paste0(substr(name, 1, 32), "..."))
      }
    }

    return(name)
  }

  all_loadings <- all_loadings %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      Feature_display = format_feature_name(Feature, View)
    ) %>%
    dplyr::ungroup()

  # Update Feature_label with formatted name (no asterisks)
  all_loadings <- all_loadings %>%
    dplyr::mutate(Feature_label = Feature_display)

  # --------------------------------------------------------------------------
  # 5. Add variance explained to factor labels
  # --------------------------------------------------------------------------
  # Calculate total variance per factor
  factor_var <- apply(var_explained, 2, sum)
  factor_labels <- setNames(
    paste0("Factor ", seq_along(factor_var), " (", round(factor_var, 1), "%)"),
    paste0("Factor", seq_along(factor_var))
  )

  all_loadings <- all_loadings %>%
    dplyr::mutate(
      Factor_label = factor_labels[Factor],
      View_label = view_labels[View]
    )

  # --------------------------------------------------------------------------
  # 6. Order features within each panel by loading value
  # --------------------------------------------------------------------------
  # Create unique ordering per View (for proper facet ordering)
  # Using tidytext-style reorder_within approach
  all_loadings <- all_loadings %>%
    dplyr::group_by(Factor, View) %>%
    dplyr::arrange(Loading, .by_group = TRUE) %>%
    dplyr::mutate(
      y_order = dplyr::row_number(),
      # Create unique feature label per view for proper ordering
      Feature_ordered = paste0(Feature_label, "___", View)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Feature_ordered = stats::reorder(Feature_ordered, Loading)
    )

  # --------------------------------------------------------------------------
  # 7. Create the diverging lollipop plot
  # --------------------------------------------------------------------------

  # Color palette (user-specified colors)
  direction_colors <- c(
    "BV-enriched" = "#F68B33",        # Orange for BV
    "Control-enriched" = "#388ECC"    # Blue for Control
  )

  # Determine facet structure based on figure type
  # Views as ROWS for independent y-axis ordering per view
  if (figure_type == "main") {
    # Views as rows, single column
    facet_formula <- View_label ~ .
    fig_width <- 6
    fig_height <- 10
  } else {
    # Views as rows, factors as columns
    facet_formula <- View_label ~ Factor_label
    fig_width <- 8
    fig_height <- 10
  }

  # Set view order
  all_loadings$View_label <- factor(
    all_loadings$View_label,
    levels = c("Taxonomy", "Metabolites (ESI+)", "Metabolites (ESI\u2212)", "Pathways")
  )

  # Helper function to clean labels (remove ___View suffix)
  clean_label <- function(x) {
    gsub("___.*$", "", x)
  }

  # Build the plot
  p <- ggplot2::ggplot(
    all_loadings,
    ggplot2::aes(
      x = Loading,
      y = Feature_ordered,
      color = Direction
    )
  ) +
    # Lollipop segments
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = Loading, yend = Feature_ordered),
      linewidth = 0.6,
      alpha = 0.8
    ) +
    # Lollipop points
    ggplot2::geom_point(size = 2.5, alpha = 0.9) +
    # Vertical reference line at 0
    ggplot2::geom_vline(
      xintercept = 0, linetype = "solid", color = "gray40", linewidth = 0.4
    ) +
    # Faceting
    ggplot2::facet_grid(
      facet_formula,
      scales = "free_y",
      space = "free_y"
    ) +
    # Colors
    ggplot2::scale_color_manual(
      values = direction_colors,
      name = "Association"
    ) +
    # X-axis with symmetric limits and 0.5 interval tick marks
    {
      # Calculate symmetric axis limits based on max absolute loading
      max_abs <- ceiling(max(abs(all_loadings$Loading), na.rm = TRUE) * 2) / 2
      ggplot2::scale_x_continuous(
        breaks = seq(-max_abs, max_abs, by = 0.5),
        limits = c(-max_abs, max_abs),
        labels = function(x) sprintf("%.1f", x)
      )
    } +
    # Clean y-axis labels (remove ___View suffix)
    ggplot2::scale_y_discrete(labels = clean_label) +
    # Labels (no title/subtitle)
    ggplot2::labs(
      x = "Loading",
      y = NULL
    ) +
    # Theme
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "gray95", color = "gray70"),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggtext::element_markdown(size = 8),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold")
    )

  # --------------------------------------------------------------------------
  # 8. Save figures
  # --------------------------------------------------------------------------
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Determine filename
  if (figure_type == "main") {
    filename_base <- "MOFA2_Factor1_loadings"
  } else {
    filename_base <- "MOFA2_loadings_supplementary"
  }

  # Save PDF
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(filename_base, ".pdf")),
    plot = p,
    width = fig_width,
    height = fig_height,
    units = "in"
  )

  # Save PNG at 300 DPI
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(filename_base, ".png")),
    plot = p,
    width = fig_width,
    height = fig_height,
    units = "in",
    dpi = 300
  )

  cat(sprintf("\nFigures saved to:\n  - %s.pdf\n  - %s.png\n",
              file.path(output_dir, filename_base),
              file.path(output_dir, filename_base)))

  # --------------------------------------------------------------------------
  # 9. Return results
  # --------------------------------------------------------------------------
  return(list(
    plot = p,
    data = all_loadings,
    factor_variance = factor_var[factors],
    diablo_overlap_count = sum(all_loadings$DIABLO_overlap)
  ))
}


# Enhanced pH analysis with metabolomics integration
analyze_enhanced_ph_patterns <- function(sample_metadata, MOFAobject) {

  cat("\n=== ENHANCED pH PATTERNS WITH MULTI-OMICS ===\n")

  # Basic pH analysis by clinical diagnosis
  ph_summary <- sample_metadata %>%
    group_by(Group) %>%
    summarise(
      Count = n(),
      Mean_pH = round(mean(Vaginal_pH, na.rm = TRUE), 2),
      SD_pH = round(sd(Vaginal_pH, na.rm = TRUE), 3),
      Min_pH = round(min(Vaginal_pH, na.rm = TRUE), 2),
      Max_pH = round(max(Vaginal_pH, na.rm = TRUE), 2),
      .groups = 'drop'
    ) %>%
    arrange(Mean_pH)

  cat("pH patterns by clinical diagnosis:\n")
  print(ph_summary)



  # MOFA factor correlation with pH
  factors <- get_factors(MOFAobject)[[1]]
  ph_values <- sample_metadata$Vaginal_pH[!is.na(sample_metadata$Vaginal_pH)]
  factor_ph_cors <- apply(factors[!is.na(sample_metadata$Vaginal_pH), ], 2,
                          function(x) cor(x, ph_values, use = "complete.obs"))

  cat("\nMOFA Factor correlations with vaginal pH:\n")
  print(round(factor_ph_cors, 3))

  # Create enhanced pH visualization
  p_ph_enhanced <- ggplot(sample_metadata, aes(x = Group, y = Vaginal_pH,
                                               fill = Group)) +
    geom_violin(alpha = 0.7, width = 0.8) +
    geom_boxplot(width = 0.3, alpha = 0.8, outlier.size = 1) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    scale_fill_manual(values = c("Control" = "#4CAF50",
                                 "Vaginosis" = "#F44336")) +
    geom_hline(yintercept = 4.5, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_hline(yintercept = 4.0, linetype = "dashed", color = "blue", alpha = 0.7) +
    labs(title = "Vaginal pH Distribution by Clinical Diagnosis",
         subtitle = "Blue line: Optimal pH (4.0), Red line: Normal threshold (4.5)",
         x = "Clinical Diagnosis",
         y = "Vaginal pH",
         fill = "Diagnosis") +
    theme_bw() +
    theme(legend.position = "none")

  ggsave("figures/MOFA2_pH_analysis.png", p_ph_enhanced)

  print(p_ph_enhanced)

  return(list(
    ph_summary = ph_summary,
    factor_ph_correlations = factor_ph_cors,
    ph_plot = p_ph_enhanced
  ))
}

# ============================================================================
# MAIN ENHANCED ANALYSIS FUNCTION
# ============================================================================

run_complete_4view_vaginal_analysis <- function(MOFAobject, sample_metadata) {

  cat("====================================================================\n")
  cat("COMPLETE 4-VIEW VAGINAL MICROBIOME ANALYSIS\n")
  cat("Metabolomics (Pos/Neg) + Pathways + Taxonomy Integration\n")
  cat("====================================================================\n")

  # Step 1: Enhanced variance analysis
  cat("\n[STEP 1] Analyzing 4-View Variance Explained...\n")
  variance_results <- analyze_4view_variance_explained(MOFAobject)

  # Step 2: Enhanced factor analysis
  cat("\n[STEP 2] Enhanced Factor Analysis...\n")
  factor1_results <- analyze_enhanced_factor_loadings(MOFAobject, factor_num = 1, top_n = 30)
  factor2_results <- analyze_enhanced_factor_loadings(MOFAobject, factor_num = 2, top_n = 20)
  factor3_results <- analyze_enhanced_factor_loadings(MOFAobject, factor_num = 3, top_n = 15)
  str(factor1_results)

  factor1_results |>
    dplyr::filter(View=="taxonomy",
                  Loading >0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  factor1_results |>
    dplyr::filter(View=="taxonomy",
                  Loading <0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  factor1_results |>
    dplyr::filter(View=="metabolomics_pos",
                  Loading >0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  factor1_results |>
    dplyr::filter(View=="metabolomics_pos",
                  Loading <0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  factor1_results |>
    dplyr::filter(View=="metabolomics_neg",
                  Loading >0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  factor1_results |>
    dplyr::filter(View=="metabolomics_neg",
                  Loading <0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  factor1_results |>
    dplyr::filter(View=="pathways",
                  Loading >0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  factor1_results |>
    dplyr::filter(View=="pathways",
                  Loading <0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  factor2_results |>
    dplyr::filter(View=="taxonomy",
                  Loading >0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  factor2_results |>
    dplyr::filter(View=="taxonomy",
                  Loading <0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  factor2_results |>
    dplyr::filter(View=="metabolomics_pos",
                  Loading >0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  factor2_results |>
    dplyr::filter(View=="metabolomics_pos",
                  Loading <0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  factor2_results |>
    dplyr::filter(View=="metabolomics_neg",
                  Loading >0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  factor2_results |>
    dplyr::filter(View=="metabolomics_neg",
                  Loading <0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  factor2_results |>
    dplyr::filter(View=="pathways",
                  Loading >0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  factor2_results |>
    dplyr::filter(View=="pathways",
                  Loading <0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  analyze_enhanced_factor_loadings(MOFAobject, factor_num = 5, top_n = 15) |>
    dplyr::filter(View=="taxonomy",
                  Loading >0) |>
    dplyr::arrange(desc(Abs_Loading)) |>
    dplyr::select(-c(Feature,View, Abs_Loading)) |>
    head(5)

  load("objects/diablo_selected_features.rda")


  glimpse(all_selected_features)


  # ---- DIABLO sectioned TXT export (3 components) -----------------------------
  # Input:  all_selected_features with cols:
  #         Feature (chr), Loading (dbl), Component (dbl/int),
  #         Block (chr: e.g., "Metab_pos","Metab_neg","Taxa","Path"),
  #         Abs_Loading (dbl)

  # Load packages (prefixed usage below)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("readr", quietly = TRUE)
  requireNamespace("stringr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  requireNamespace("purrr", quietly = TRUE)

  # ---- Config ------------------------------------------------------------------
  top_n  <- 5

  # Map mixOmics DIABLO block names -> section labels used in the report
  block_map <- c(
    "Metab_pos" = "metabolomics_pos",
    "Metab_neg" = "metabolomics_neg",
    "Taxa"      = "taxonomy",
    "Path"      = "pathways",
    "Pathway"   = "pathways",
    "Pathways"  = "pathways"
  )

  # Preferred display order of sections (if present)
  view_order <- c("taxonomy", "metabolomics_pos", "metabolomics_neg", "pathways")

  # ---- Prep --------------------------------------------------------------------
  df <- all_selected_features %>%
    dplyr::mutate(
      View = dplyr::recode(Block, !!!block_map),
      Abs_Loading = if (!"Abs_Loading" %in% names(all_selected_features))
        abs(Loading) else Abs_Loading
    ) %>%
    dplyr::filter(!is.na(View)) %>%
    dplyr::mutate(Component = as.integer(Component))

  # Helper to format one (view x direction) into lines
  format_section_table <- function(data, section_title) {
    # data must contain Feature, Loading
    header <- c(section_title, "Feature\tLoading")
    if (nrow(data) == 0L) {
      return(c(header, "(no features)\n"))
    }
    rows <- data %>%
      dplyr::mutate(
        Loading = sprintf("%.6f", Loading)
      ) %>%
      dplyr::transmute(line = paste0(Feature, "\t", Loading)) %>%
      dplyr::pull(line)
    c(header, rows, "")  # blank line after table
  }

  # ---- Build the full report ---------------------------------------------------
  all_lines <- character(0)

  for (k in sort(unique(df$Component))) {
    all_lines <- c(all_lines, paste0("Component ", k))

    comp_df <- df %>%
      dplyr::filter(Component == k)

    # iterate views in preferred order, but only those present
    views_here <- comp_df %>%
      dplyr::distinct(View) %>%
      dplyr::pull(View)
    views_here <- intersect(view_order, views_here)

    for (vv in views_here) {
      sub <- comp_df %>% dplyr::filter(View == vv)

      # Enriched (Loading > 0)
      tab_pos <- sub %>%
        dplyr::filter(Loading > 0) %>%
        dplyr::arrange(dplyr::desc(Abs_Loading)) %>%
        dplyr::select(Feature, Loading) %>%
        dplyr::slice_head(n = top_n)

      all_lines <- c(
        all_lines,
        format_section_table(tab_pos, paste0("DIABLO - ", vv, " - enriched"))
      )

      # Depleted (Loading < 0)
      tab_neg <- sub %>%
        dplyr::filter(Loading < 0) %>%
        dplyr::arrange(dplyr::desc(Abs_Loading)) %>%
        dplyr::select(Feature, Loading) %>%
        dplyr::slice_head(n = top_n)

      all_lines <- c(
        all_lines,
        format_section_table(tab_neg, paste0("DIABLO - ", vv, " - depleted"))
      )
    }

    # spacer between components
    all_lines <- c(all_lines, "")
  }





  # gv <- data.table::fread("data/gv_mags.tsv")
  #
  # gv |>
  #   dplyr::mutate(
  #     spp = stringr::str_extract(mags, "(?<=_)[^_]+.*"),     # everything after the first "_"
  #     Sample = stringr::str_extract(mags, "^[^_]+")           # everything before the first "_"
  #   ) |>
  #   dplyr::arrange(desc(mags)) |>
  #   dplyr::select(Sample, spp) |>
  #   group_by(spp) |>
  #   summarise(count= n())



  # Step 3: Enhanced microbiome state classification
  cat("\n[STEP 3] Enhanced Microbiome State Classification...\n")
  sample_metadata <- classify_enhanced_microbiome_states(MOFAobject, sample_metadata)

  # Step 4: Enhanced classification analysis
  cat("\n[STEP 4] Enhanced 4-View Classification...\n")
  classification_results <- run_enhanced_4view_classification(MOFAobject, sample_metadata)

  # Step 5: Enhanced pH analysis
  cat("\n[STEP 5] Enhanced pH Analysis...\n")
  ph_results <- analyze_enhanced_ph_patterns(sample_metadata, MOFAobject)

  # Step 6: Enhanced visualizations
  cat("\n[STEP 6] Creating Enhanced Visualizations...\n")
  heatmap_results <- create_enhanced_4view_heatmaps(MOFAobject, top_n = 20)

  # Step 7: Clinical concordance analysis
  cat("\n[STEP 7] Clinical-Microbiome Concordance...\n")
  cross_tab <- table(sample_metadata$Group, sample_metadata$Enhanced_Microbiome_State)
  cat("Clinical vs Enhanced Microbiome State Cross-tabulation:\n")
  print(cross_tab)

  # Calculate key metrics
  control_total <- sum(sample_metadata$Group == "Control")
  control_optimal <- sum(sample_metadata$Group == "Control" &
                           sample_metadata$Enhanced_Microbiome_State == "Optimal_Lactobacilli")
  control_bv <- sum(sample_metadata$Group == "Control" &
                      sample_metadata$Enhanced_Microbiome_State %in% c("Moderate_BV", "Severe_BV"))

  cat("\n=== KEY DISCOVERIES ===\n")
  cat(sprintf("1. OPTIMAL HEALTH: %d/%d controls (%.1f%%) have optimal microbiome\n",
              control_optimal, control_total, control_optimal/control_total*100))
  cat(sprintf("2. ASYMPTOMATIC BV: %d/%d controls (%.1f%%) have BV biology\n",
              control_bv, control_total, control_bv/control_total*100))

  # Return comprehensive results
  return(list(
    variance_results = variance_results,
    factor1_results = factor1_results,
    factor2_results = factor2_results,
    factor3_results = factor3_results,
    sample_metadata = sample_metadata,
    classification_results = classification_results,
    ph_results = ph_results,
    heatmap_results = heatmap_results,
    cross_tab = cross_tab
  ))
}

# ============================================================================
# EXECUTION
# ============================================================================

cat("\n====================================================================\n")
cat("READY TO RUN COMPLETE 4-VIEW ANALYSIS\n")
cat("Execute: analysis_results <- run_complete_4view_vaginal_analysis(MOFAobject, sample_metadata)\n")
cat("====================================================================\n")

# USAGE:
analysis_results <- run_complete_4view_vaginal_analysis(MOFAobject, sample_metadata)

# 1) Which factors align with BV and vaginal pH?
Z <- MOFA2::get_factors(MOFAobject, as.data.frame = FALSE)$group1  # samples x factors
W_tax <- MOFA2::get_weights(MOFAobject, views = "taxonomy", as.data.frame = FALSE)[[1]] # taxa x factors

# Correlate factors with covariates (BV coded 0/1; add pH, etc.)
covar_df <- data.frame(
  sample = rownames(Z),
  BV = as.numeric(Y_aligned[rownames(Z)] == "BV"),
  Vaginal_pH = sample_metadata_aligned[rownames(Z), "Vaginal_pH"]
)
rownames(covar_df) <- covar_df$sample; covar_df$sample <- NULL

MOFA2::correlate_factors_with_covariates(
  MOFAobject, covariates = covar_df, factors = "all", plot = "r"
)

# 2) Flip factor signs so that positive Z means "more BV"
flip <- sign(stats::cor(Z, covar_df$BV, use = "pairwise.complete.obs"))  # length K
flip[is.na(flip)] <- 1
Z_aligned    <- sweep(Z, 2, flip, `*`)
W_tax_aligned <- sweep(W_tax, 2, flip, `*`)

# 3) Per-category enrichment of taxa on each factor (signed weights)
library(dplyr); library(tidyr); library(tibble)

# Auto-categorize taxa based on species names
categorize_taxon <- function(taxon_name) {
  name_lower <- tolower(taxon_name)
  dplyr::case_when(
    grepl("lactobacillus.*iners|l\\..*iners", name_lower) ~ "Transitional_Lactobacilli",
    grepl("lactobacillus|limosilactobacillus|lacticaseibacillus", name_lower) ~ "Protective_Lactobacilli",
    grepl("gardnerella|fannyhessea|atopobium", name_lower) ~ "Classic_BV_Pathogens",
    grepl("prevotella|sneathia|dialister|megasphaera|mobiluncus", name_lower) ~ "BV_Associated_Anaerobes",
    grepl("veillonella|coriobacter", name_lower) ~ "BV_Associated_Anaerobes",
    grepl("bifidobacterium|corynebacterium", name_lower) ~ "Commensal_Bacteria",
    grepl("streptococcus|staphylococcus|enterococcus", name_lower) ~ "Opportunistic_Pathogens",
    grepl("escherichia|klebsiella|enterobacter", name_lower) ~ "Enterobacteriaceae",
    grepl("candida|saccharomyces", name_lower) ~ "Fungi",
    TRUE ~ "Other_Taxa"
  )
}

taxa_cat <- tibble::tibble(
  taxon = rownames(W_tax_aligned),
  Category = sapply(rownames(W_tax_aligned), categorize_taxon)
)

enrich <- lapply(seq_len(ncol(W_tax_aligned)), function(k){
  dplyr::tibble(taxon = rownames(W_tax_aligned),
                weight = W_tax_aligned[, k]) |>
    dplyr::left_join(taxa_cat, by = "taxon") |>
    dplyr::group_by(Category) |>
    dplyr::summarise(sum_weight = sum(weight), .groups = "drop") |>
    dplyr::mutate(Factor = paste0("F", k))
}) |> dplyr::bind_rows()

print(enrich)  # positive sum_weight -> category increases with BV (after sign alignment)

# 4) How much of taxonomy variance each factor explains
r2 <- MOFA2::calculate_variance_explained(MOFAobject)
MOFA2::plot_variance_explained(MOFAobject, x = "view", y = "factor", factors = 1:7)

# ============================================================================
# GENERATE MOFA LOADINGS DIVERGING LOLLIPOP FIGURES
# ============================================================================

# Extract DIABLO features for overlap annotation
load("objects/diablo_selected_features.rda")
diablo_features_all <- unique(all_selected_features$Feature)
cat(sprintf("\nFound %d unique DIABLO features for overlap annotation\n",
            length(diablo_features_all)))

# Generate MAIN TEXT figure: Factor 1 only
cat("\n--- Generating Main Text Figure (Factor 1) ---\n")
mofa_fig_main <- create_mofa_loadings_figure(
  MOFAobject = MOFAobject,
  factors = c(1),
  top_n = 5,
  diablo_features = diablo_features_all,
  figure_type = "main",
  output_dir = "figures"
)
cat(sprintf("DIABLO overlap in main figure: %d features\n",
            mofa_fig_main$diablo_overlap_count))

# Generate SUPPLEMENTARY figure: Factor 1 + Factor 2
cat("\n--- Generating Supplementary Figure (Factors 1 + 2) ---\n")
mofa_fig_supp <- create_mofa_loadings_figure(
  MOFAobject = MOFAobject,
  factors = c(1, 2),
  top_n = 5,
  diablo_features = diablo_features_all,
  figure_type = "supplementary",
  output_dir = "figures"
)
cat(sprintf("DIABLO overlap in supplementary figure: %d features\n",
            mofa_fig_supp$diablo_overlap_count))

# Print the main figure
print(mofa_fig_main$plot)

cat("\n=== MOFA LOADINGS FIGURES COMPLETE ===\n")
cat("Main figure: figures/MOFA2_Factor1_loadings.pdf/png\n")
cat("Supplementary: figures/MOFA2_loadings_supplementary.pdf/png\n")

# ============================================================================
# CREATE SUPPLEMENTARY TABLES S9 AND S10 (Excel workbooks)
# ============================================================================

cat("\n=== CREATING SUPPLEMENTARY TABLES (Excel) ===\n")

# --------------------------------------------------------------------------
# Table S9: Variance Explained by MOFA2 Factors
# --------------------------------------------------------------------------
var_explained <- MOFA2::get_variance_explained(MOFAobject)$r2_per_factor$group1

# Convert to tidy format for Table S9
table_s9 <- as.data.frame(var_explained) %>%
  tibble::rownames_to_column("View") %>%
  dplyr::mutate(
    View = dplyr::case_when(
      View == "metabolomics_pos" ~ "Metabolomics (ESI+)",
      View == "metabolomics_neg" ~ "Metabolomics (ESI\u2212)",
      View == "pathways" ~ "Functional Pathways",
      View == "taxonomy" ~ "Taxonomy",
      TRUE ~ View
    )
  )

# Add total variance per factor (sum across views)
total_row <- c("Total (all views)", colSums(var_explained))
names(total_row) <- names(table_s9)
table_s9 <- rbind(table_s9, total_row)

# Round numeric columns
numeric_cols <- setdiff(names(table_s9), "View")
table_s9[numeric_cols] <- lapply(table_s9[numeric_cols], function(x) {
  round(as.numeric(x), 2)
})

cat("Table S9 - Variance explained by MOFA2 factors:\n")
print(table_s9)

# --------------------------------------------------------------------------
# Table S10: Complete Feature Loadings for All Factors
# --------------------------------------------------------------------------
weights <- MOFA2::get_weights(MOFAobject)

# Build comprehensive loadings table
table_s10 <- purrr::map_dfr(names(weights), function(view_name) {
  view_weights <- weights[[view_name]]

  # Convert to data frame
  df <- as.data.frame(view_weights) %>%
    tibble::rownames_to_column("Feature") %>%
    dplyr::mutate(
      View = dplyr::case_when(
        view_name == "metabolomics_pos" ~ "Metabolomics (ESI+)",
        view_name == "metabolomics_neg" ~ "Metabolomics (ESI\u2212)",
        view_name == "pathways" ~ "Functional Pathways",
        view_name == "taxonomy" ~ "Taxonomy",
        TRUE ~ view_name
      )
    ) %>%
    dplyr::select(View, Feature, dplyr::everything())

  return(df)
})

# Rename factor columns for clarity
factor_cols <- grep("^Factor", names(table_s10), value = TRUE)
for (fc in factor_cols) {
  table_s10[[fc]] <- round(table_s10[[fc]], 4)
}

# Add max absolute loading column for sorting/filtering
table_s10 <- table_s10 %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Max_Abs_Loading = max(abs(dplyr::c_across(dplyr::starts_with("Factor"))))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(View, dplyr::desc(Max_Abs_Loading))

# Add DIABLO overlap indicator
table_s10 <- table_s10 %>%
  dplyr::mutate(
    DIABLO_Overlap = dplyr::if_else(Feature %in% diablo_features_all, "Yes", "No")
  )

cat(sprintf("\nTable S10 - Complete feature loadings: %d features\n", nrow(table_s10)))
cat(sprintf("  DIABLO overlap features: %d\n",
            sum(table_s10$DIABLO_Overlap == "Yes")))

# --------------------------------------------------------------------------
# Create Separate Excel Workbooks (one per table)
# --------------------------------------------------------------------------
header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  border = "bottom",
  halign = "center"
)

# --- Table S9 Workbook (Variance Explained) ---
wb_s9 <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_s9, "Table_S9_Variance")

# Add title
openxlsx::writeData(wb_s9, "Table_S9_Variance",
                    "Table S9. Variance explained (%) by MOFA2 factors across omics views",
                    startRow = 1, startCol = 1)

# Add data
openxlsx::writeData(wb_s9, "Table_S9_Variance", table_s9, startRow = 3)

# Style header
openxlsx::addStyle(wb_s9, "Table_S9_Variance", header_style,
                   rows = 3, cols = 1:ncol(table_s9), gridExpand = TRUE)

# Auto-width columns
openxlsx::setColWidths(wb_s9, "Table_S9_Variance",
                       cols = 1:ncol(table_s9), widths = "auto")

output_s9 <- "supplementary_tables/Table_S9.xlsx"
openxlsx::saveWorkbook(wb_s9, output_s9, overwrite = TRUE)

# --- Table S10 Workbook (Complete Loadings) ---
wb_s10 <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_s10, "Table_S10_Loadings")

# Add title
openxlsx::writeData(wb_s10, "Table_S10_Loadings",
                    "Table S10. Complete MOFA2 factor loadings for all features",
                    startRow = 1, startCol = 1)

# Add data
openxlsx::writeData(wb_s10, "Table_S10_Loadings", table_s10, startRow = 3)

# Style header
openxlsx::addStyle(wb_s10, "Table_S10_Loadings", header_style,
                   rows = 3, cols = 1:ncol(table_s10), gridExpand = TRUE)

# Conditional formatting for DIABLO overlap
diablo_style <- openxlsx::createStyle(fgFill = "#FFFACD")  # Light yellow
openxlsx::conditionalFormatting(
  wb_s10, "Table_S10_Loadings",
  cols = 1:ncol(table_s10),
  rows = 4:(nrow(table_s10) + 3),
  rule = '$N4="Yes"',
  type = "expression",
  style = diablo_style
)

# Auto-width columns
openxlsx::setColWidths(wb_s10, "Table_S10_Loadings",
                       cols = 1:ncol(table_s10), widths = "auto")

output_s10 <- "supplementary_tables/Table_S10.xlsx"
openxlsx::saveWorkbook(wb_s10, output_s10, overwrite = TRUE)

cat(sprintf("\n=== SUPPLEMENTARY TABLES SAVED ===\n"))
cat(sprintf("Table S9: %s\n", output_s9))
cat(sprintf("Table S10: %s\n", output_s10))
