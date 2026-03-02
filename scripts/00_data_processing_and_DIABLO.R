# ============================================================================
# 00_data_processing_and_DIABLO.R
# Title: MitCH Study: Data Processing, Differential Abundance, and DIABLO Integration
#
# Purpose:
#   Processes raw sequencing + metabolomics data into multi-omics objects for
#   integration. This script performs:
#     - Read count and MetaPhlAn4 species profile merging across flowcells
#     - Phyloseq object construction with decontam filtering
#     - ANCOM-BC2 taxonomic differential abundance (BV vs Control)
#     - Functional pathway differential abundance (ANCOM-BC2)
#     - Alpha and beta diversity analysis (Shannon, Observed, Evenness; Bray-Curtis)
#     - LC-MS metabolomics normalization (PQN, kNN imputation, batch correction)
#     - Limma differential analysis for positive and negative ionization modes
#     - CLR transformation and sample alignment for multi-omics integration
#     - DIABLO design matrix specification, model fitting, and visualization
#     - Gardnerella pangenome analysis with gene-metabolite correlations
#
# External data dependencies:
#   Raw reads deposited at ENA under accession PRJEB108308 (available upon
#   publication). Metabolomics and clinical metadata available from the
#   corresponding authors upon reasonable request.
#
# Notes:
#   - Downstream scripts (01, 02) run independently using cached objects
#     in objects/ and metadata/
#   - This script requires raw data files not included in this repository
#
# ============================================================================
#
# DATA_SOURCES
# ============================================================================
# The following external files are required to run this script. They are NOT
# included in the repository due to data-sharing restrictions or because they
# are generated from upstream pipelines.
#
# INPUT FILES:
#
#   data/MatchMetaData.SeqMetaDATA.xlsx
#     -> Clinical metadata, available from authors upon request
#
#   data/E200019361_read_counts.txt
#   data/E250045577_read_counts.txt
#   data/E250045590_read_counts.txt
#     -> Read count data from StaG pipeline, generated from ENA raw reads (PRJEB108308)
#
#   data/E200019361_species.tsv
#   data/E250045577_species.tsv
#   data/E250045590_species.tsv
#     -> MetaPhlAn4 species profiles from StaG pipeline
#
#   data/polar_negative_results.tsv
#   data/polar_positive_results.tsv
#     -> LC-MS metabolomics data, available from authors upon request
#
#   (Gardnerella pangenome data used in scripts/02_genome_metabolome.R)
#
#   objects/merged_path_abundance.rda
#     -> HUMAnN3 pathway abundance profiles (generated from StaG pipeline)
#
#   objects/merged_pathway_coverage.rda
#     -> HUMAnN3 pathway coverage profiles (generated from StaG pipeline)
#
#   objects/merged_gene_family.rda
#     -> HUMAnN3 gene family profiles (generated from StaG pipeline)
#
#   objects/taxa_ancombc2_adjust_full.rda
#   objects/taxa_ancombc2_basic.rda
#     -> Pre-computed ANCOM-BC2 taxonomic DA results (adjusted and basic models)
#
#   objects/path_ancombc2_basic.rda
#   objects/path_ancombc2_adjusted.rda
#     -> Pre-computed ANCOM-BC2 pathway DA results
#
#   objects/tune_ncomp_full_design_BV_v3.rda
#   objects/tune_ncomp_biological_design_BV_v3.rda
#   objects/tune_ncomp_hybrid_design_BV_v3.rda
#   objects/tune_ncomp_datadriven_design_BV_v3.rda
#   objects/tune_ncomp_null_design_BV_v3.rda
#     -> Pre-computed DIABLO cross-validation tuning results
#
#
# OUTPUT FILES (saved to repository):
#
#   objects/multiomics_views.rda          (multi-omics data list: X)
#   objects/group_labels.rda              (outcome vector: Y_aligned)
#   metadata/sample_metadata.rda          (aligned metadata: meta_data)
#   objects/bact_object.rda               (phyloseq-ready bacteriome objects)
#   objects/BACT_complete_meta.rda        (complete metadata with sample types)
#   objects/bv_path_res.rda               (pathway DA results)
#   figures/volcano_da_bv_tax_BV.png
#   figures/volcano_da_bv_tax_BV_total.png
#   figures/volcano_da_path_BV_full.png
#   figures/alpha_diversity_panel_plot.png
#   figures/Figure1.png
#   figures/volcano_da_bv_pos_BV_adj_v2.png
#   figures/volcano_da_bv_neg_BV__adj_v2.png
#   figures/combined_volcano_bv_pos_neg_v2.png
#   figures/integration_plots/DIABLO_sample_plot_comp1-2.png
#   figures/integration_plots/DIABLO_AUROC_comprehensive.pdf
#   figures/integration_plots/circos_plot_diablo_v4.pdf
#   figures/integration_plots/circos_plot_diablo_v3.png
#   supplementary_tables/Table_S2.xlsx    (PERMANOVA results)
#   supplementary_tables/Table_S3.xlsx    (taxonomic DA, ANCOM-BC2)
#   supplementary_tables/Table_S4.xlsx    (pathway DA, ANCOM-BC2)
#   supplementary_tables/Table_S5.xlsx    (metabolomics positive mode, limma)
#   supplementary_tables/Table_S6.xlsx    (metabolomics negative mode, limma)
#   supplementary_tables/Table_S7.xlsx    (DIABLO variance explained)
#   supplementary_tables/Table_S8.xlsx    (DIABLO feature loadings)
#   (Table_S11 generated by scripts/02_genome_metabolome.R)
#
# ============================================================================

# --- 1. SETUP AND LIBRARY LOADING -------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(readxl)
  library(reshape2)
  library(vegan)
  library(ggplot2)
  library(compositions)
  library(stringr)
  library(scales)
  library(pheatmap)
  library(ggrepel)
  library(DESeq2)
  library(glue)
  library(ANCOMBC)
  library(decontam)
  library(phyloseq)
  library(limma)
  library(randomForest)
  library(caret)
  library(mixOmics)
  library(tibble)
  library(magrittr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(matrixStats)
  library(impute)
  library(sva)
  library(patchwork)
  library(openxlsx)
  library(data.table)
  library(RColorBrewer)
  library(broom)
  library(BiocParallel)
  library(parallel)
  library(gridExtra)
  library(viridis)
  library(ade4)
})

sessionInfo()

# --- 2. DATA LOADING AND PREPROCESSING --------------------------------------

# 2.1 Metadata loading
# NOTE: requires raw data not included in repository
MetaData <- readxl::read_excel("data/MatchMetaData.SeqMetaDATA.xlsx")

# 2.2 Read count processing across flowcells
process_read_data <- function(flowcell_id) {
  # NOTE: requires raw data not included in repository
  file_path <- paste0("data/", flowcell_id, "_read_counts.txt")
  read_data <- read.table(file_path, sep = "\t", header = TRUE)

  processed_data <- read_data %>%
    dplyr::mutate(Sample = gsub("236_Mitch__", "", Sample)) %>%
    dplyr::mutate(Sample = gsub(paste0("__", flowcell_id, "_L1.*"), "", Sample)) %>%
    dplyr::mutate(Sample = dplyr::case_when(
      grepl("^neg_ex", Sample) ~ paste0("neg_ex__", flowcell_id),
      grepl("^neg_seq", Sample) ~ paste0("neg_seq__", flowcell_id),
      grepl("^zymomock_ex", Sample) ~ paste0("zymomock_ex__", flowcell_id),
      grepl("^zymomock_seq", Sample) ~ paste0("zymomock_seq__", flowcell_id),
      TRUE ~ Sample
    )) %>%
    dplyr::mutate(reads = after_kraken2_host_removal) %>%
    dplyr::select(Sample, reads)

  processed_data <- cbind(processed_data, flowcell_id)
  return(processed_data)
}

flowcell_id <- c("E200019361", "E250045577", "E250045590")
processed_read_data_list <- lapply(flowcell_id, process_read_data)
read_data <- do.call(rbind, processed_read_data_list)

# Merge with metadata
merged_read_data <- read_data %>%
  dplyr::filter(Sample %in% MetaData$SeqID) %>%
  dplyr::left_join(MetaData, by = c("Sample" = "SeqID")) %>%
  dplyr::select(-flowcell_id)

# 2.3 MetaPhlAn species-level data processing
process_metaphlan_spp_data <- function(flowcell_id) {
  # NOTE: requires raw data not included in repository
  file_path <- paste0("data/", flowcell_id, "_species.tsv")
  data <- read.table(file_path, sep = "\t", header = TRUE)

  processed_data <- data %>%
    dplyr::filter(grepl("s__", clade_name)) %>%
    dplyr::filter(!grepl("t__", clade_name)) %>%
    dplyr::mutate(clade_name = gsub(".*s__", "", clade_name))

  processed_data <- tibble::column_to_rownames(processed_data, var = "clade_name")

  regex <- paste0("(F[A-Z][0-9]+|neg_ex__", flowcell_id, "|neg_seq__", flowcell_id,
                  "|zymomock_ex__", flowcell_id, "|zymomock_seq__", flowcell_id, ")")

  colnames(processed_data) <- stringr::str_extract(colnames(processed_data), regex)
  processed_data <- tibble::rownames_to_column(processed_data, var = "clade_name")
  return(processed_data)
}

processed_spp_data_list <- lapply(flowcell_id, process_metaphlan_spp_data)
merged_spp_data <- Reduce(function(x, y) dplyr::left_join(x, y, by = "clade_name"), processed_spp_data_list)
merged_spp_data[is.na(merged_spp_data)] <- 0

# 2.4 Pseudo-count conversion from relative abundance
read_selected <- read_data %>%
  dplyr::select(Sample, reads)

merged_spp_counts <- merged_spp_data
read_selected <- read_selected[match(colnames(merged_spp_data)[-1], read_selected$Sample), ]
merged_spp_counts[, -1] <- round(sweep(merged_spp_data[, -1], 2, read_selected$reads / 100, "*"))

# Load pre-processed metadata (public version included in repository)
load("metadata/sample_metadata.rda")
Bact_meta <- meta_data

# Derive BACT_complete_meta from public metadata
BACT_complete_meta <- Bact_meta %>%
  dplyr::select(Sample, Group, Vaginal_pH, Age, Antibiotic, BV_History)

save(BACT_complete_meta, file = "objects/BACT_complete_meta.rda")

# --- 3. PHYLOSEQ OBJECT CONSTRUCTION ----------------------------------------

# 3.1 Prepare data for phyloseq
neg_mock_cols <- grep("neg_ex|neg_seq|zymomock_ex|zymomock_seq", colnames(merged_spp_data), value = TRUE)
neg_mock_cols <- neg_mock_cols[!is.na(neg_mock_cols)]

Bact_taxaRelAbundance <- merged_spp_data %>%
  dplyr::select(clade_name, dplyr::any_of(BACT_complete_meta$Sample), dplyr::all_of(neg_mock_cols))

neg_mock_cols <- grep("neg_ex|neg_seq|zymomock_ex|zymomock_seq", colnames(merged_spp_counts), value = TRUE)
neg_mock_cols <- neg_mock_cols[!is.na(neg_mock_cols)]

Bact_taxaPseudoCounts <- merged_spp_counts %>%
  dplyr::select(clade_name, dplyr::any_of(BACT_complete_meta$Sample), dplyr::all_of(neg_mock_cols))

Bact_lib_size <- read_selected %>%
  dplyr::filter(Sample %in% colnames(Bact_taxaPseudoCounts)[-1])

# Add neg/mock samples to metadata
neg_mock_to_add <- setdiff(neg_mock_cols, Bact_meta$Sample)
if (length(neg_mock_to_add) > 0) {
  neg_mock_df <- as.data.frame(matrix(NA, nrow = length(neg_mock_to_add), ncol = ncol(Bact_meta)))
  colnames(neg_mock_df) <- colnames(Bact_meta)
  neg_mock_df$Sample <- neg_mock_to_add
  Bact_meta <- rbind(Bact_meta, neg_mock_df)
  rownames(Bact_meta) <- Bact_meta$Sample
}

# Save objects
bact_object <- list(Bact_lib_size, Bact_taxaPseudoCounts, Bact_taxaRelAbundance, Bact_meta)
names(bact_object) <- c("Bact_lib_size", "Bact_taxaPseudoCounts", "Bact_taxaRelAbundance", "Bact_meta")
save(bact_object, file = "objects/bact_object.rda")

# Update metadata with sample type
rownames(Bact_meta) <- Bact_meta$Sample
BACT_complete_meta <- Bact_meta %>%
  dplyr::mutate(
    Sample_ID = Sample,
    type = ifelse(grepl("neg_ex|neg_seq|zymomock_ex|zymomock_seq", rownames(Bact_meta)), "neg", "sample")
  ) %>%
  as.data.frame()
rownames(BACT_complete_meta) <- BACT_complete_meta$Sample

# 3.2 Create phyloseq object
stopifnot(all(colnames(Bact_taxaPseudoCounts)[-1] == BACT_complete_meta$Sample))

taxa_names <- Bact_taxaPseudoCounts$clade_name
otu_mat <- as.matrix(Bact_taxaPseudoCounts[, -1])
rownames(otu_mat) <- taxa_names

otu_physeq <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)

tax_df <- data.frame(
  Species = taxa_names,
  row.names = taxa_names,
  stringsAsFactors = FALSE
)
tax_physeq <- phyloseq::tax_table(as.matrix(tax_df))

sample_physeq <- phyloseq::sample_data(BACT_complete_meta)

physeq_bact <- phyloseq::phyloseq(otu_physeq, tax_physeq, sample_physeq)

# Set factor levels
phyloseq::sample_data(physeq_bact)$Group <- factor(phyloseq::sample_data(physeq_bact)$Group)
phyloseq::sample_data(physeq_bact)$BV_History <- factor(phyloseq::sample_data(physeq_bact)$BV_History)
phyloseq::sample_data(physeq_bact)$Antibiotic <- factor(phyloseq::sample_data(physeq_bact)$Antibiotic)
phyloseq::sample_data(physeq_bact)$type <- factor(phyloseq::sample_data(physeq_bact)$type)

phyloseq::sample_data(physeq_bact)$Group <- factor(
  phyloseq::sample_data(physeq_bact)$Group,
  levels = c("Control", "BV")
)

# 3.3 Decontam analysis
df <- as(phyloseq::sample_data(physeq_bact), "data.frame")
df$LibrarySize <- as.numeric(phyloseq::sample_sums(physeq_bact))
df <- df[order(df$LibrarySize), , drop = FALSE]
df <- df %>% dplyr::mutate(Index = seq(nrow(df)))

ggplot2::ggplot(data = df, ggplot2::aes(x = Index, y = LibrarySize, color = type)) +
  ggplot2::geom_point()

phyloseq::sample_data(physeq_bact)$is.neg <- phyloseq::sample_data(physeq_bact)$type == "neg"
contamdf.prev <- decontam::isContaminant(physeq_bact, method = "prevalence", neg = "is.neg", threshold = 0.1)
contaminant_species <- rownames(contamdf.prev)[contamdf.prev$contaminant]

taxa_counts <- phyloseq::taxa_sums(physeq_bact)
threshold_stat <- quantile(taxa_counts, probs = 0.1)
low_count_species_stat <- names(taxa_counts[taxa_counts < threshold_stat])
threshold_arbitrary <- 500
low_count_species_arbitrary <- names(taxa_counts[taxa_counts < threshold_arbitrary])

# --- HELPER FUNCTIONS -------------------------------------------------------

ensure_taxa_rows <- function(ps) {
  if (!phyloseq::taxa_are_rows(ps)) {
    ps <- phyloseq::otu_table(t(phyloseq::otu_table(ps)), taxa_are_rows = TRUE) %>%
      (\(otu) phyloseq::phyloseq(otu, phyloseq::sample_data(ps), phyloseq::tax_table(ps)))()
  }
  ps
}

prune_samples_by_pattern <- function(ps,
                                     pattern = "neg_ex|neg_seq|zymomock_ex|zymomock_seq",
                                     cols_to_check = c("Sample", "Sample_ID"),
                                     ignore_case = TRUE) {
  stopifnot(inherits(ps, "phyloseq"))
  sdf <- as.data.frame(phyloseq::sample_data(ps))

  rm_flag <- grepl(pattern, rownames(sdf), ignore.case = ignore_case)

  for (cn in cols_to_check) {
    if (cn %in% colnames(sdf)) {
      rm_flag <- rm_flag | grepl(pattern, as.character(sdf[[cn]]), ignore.case = ignore_case)
    }
  }

  n_remove <- sum(rm_flag)
  if (n_remove > 0) {
    ps <- phyloseq::prune_samples(!rm_flag, ps)
  }
  message(sprintf("Removed %d samples by pattern '%s'. Remaining: %d",
                  n_remove, pattern, phyloseq::nsamples(ps)))
  ps
}

filter_phyloseq_taxa <- function(ps,
                                 contaminant_species = character(0),
                                 min_total_reads = 500,
                                 rel_abund_thresh = 0.00005,
                                 min_samples = 10) {

  ps <- ensure_taxa_rows(ps)
  n_taxa_start <- phyloseq::ntaxa(ps)

  taxa_present <- phyloseq::taxa_names(ps)
  to_drop_contam <- intersect(taxa_present, contaminant_species)
  n_contam <- length(to_drop_contam)
  if (n_contam > 0) {
    ps <- phyloseq::prune_taxa(!(phyloseq::taxa_names(ps) %in% to_drop_contam), ps)
  }

  taxa_sums_vec <- phyloseq::taxa_sums(ps)
  low_read_taxa <- names(taxa_sums_vec[taxa_sums_vec < min_total_reads])
  if (length(low_read_taxa) > 0) {
    otu_mat <- as.matrix(phyloseq::otu_table(ps))
    rows_to_zero <- match(low_read_taxa, rownames(otu_mat))
    otu_mat[rows_to_zero, ] <- 0L
    ps <- phyloseq::phyloseq(
      phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE),
      phyloseq::sample_data(ps),
      phyloseq::tax_table(ps)
    )
  }

  zero_all <- phyloseq::taxa_sums(ps) == 0
  n_zero_all <- sum(zero_all)
  if (n_zero_all > 0) {
    ps <- phyloseq::prune_taxa(!zero_all, ps)
  }

  ps_rel <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
  rel_mat <- as.matrix(phyloseq::otu_table(ps_rel))
  keep_by_ra <- rowSums(rel_mat > rel_abund_thresh, na.rm = TRUE) >= min_samples
  n_drop_ra <- sum(!keep_by_ra)
  if (n_drop_ra > 0) {
    ps <- phyloseq::prune_taxa(keep_by_ra, ps)
  }

  n_taxa_end <- phyloseq::ntaxa(ps)
  message(
    sprintf("Start taxa: %d | Removed contaminants: %d | Zeroed-all (< %d reads): %d | ",
            n_taxa_start, n_contam, min_total_reads, n_zero_all),
    sprintf("Dropped by rel. abundance (>%.5f in >=%d samples): %d | Final taxa: %d",
            rel_abund_thresh, min_samples, n_drop_ra, n_taxa_end)
  )
  ps
}

filter_phyloseq_path <- function(ps,
                                 rel_abund_thresh = 0.00005,
                                 min_samples = 10) {

  ps <- ensure_taxa_rows(ps)
  n_taxa_start <- phyloseq::ntaxa(ps)

  ps_rel <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
  rel_mat <- as.matrix(phyloseq::otu_table(ps_rel))
  keep_by_ra <- rowSums(rel_mat > rel_abund_thresh, na.rm = TRUE) >= min_samples
  n_drop_ra <- sum(!keep_by_ra)
  if (n_drop_ra > 0) {
    ps <- phyloseq::prune_taxa(keep_by_ra, ps)
  }

  n_taxa_end <- phyloseq::ntaxa(ps)
  message(
    sprintf("Dropped by rel. abundance (>%.5f in >=%d samples): %d | Final taxa: %d",
            rel_abund_thresh, min_samples, n_drop_ra, n_taxa_end)
  )
  ps
}

# Apply filtering
physeq_bact_clean <- prune_samples_by_pattern(
  physeq_bact,
  pattern = "neg_ex|neg_seq|zymomock_ex|zymomock_seq",
  cols_to_check = c("Sample", "Sample_ID"),
  ignore_case = TRUE
)

physeq_bact_filtered <- filter_phyloseq_taxa(
  ps = physeq_bact_clean,
  contaminant_species = contaminant_species,
  min_total_reads = threshold_arbitrary,
  rel_abund_thresh = 0.00005,
  min_samples = 10
)

# --- 4. TAXONOMIC DIFFERENTIAL ABUNDANCE (ANCOM-BC2) ------------------------

# 4.1 Prepare data
sdf <- as.data.frame(phyloseq::sample_data(physeq_bact_filtered))
stopifnot("Group" %in% names(sdf))
phyloseq::sample_data(physeq_bact_filtered)[["Group"]] <-
  droplevels(factor(sdf[["Group"]], levels = c("Control", "BV")))

# 4.2 Load pre-computed ANCOM-BC2 results (full model with Age)
load("objects/taxa_ancombc2_adjust_full.rda")
load("objects/taxa_ancombc2_basic.rda")

res_raw <- taxa_ancombc2_adjust_full$res
res_raw_total <- taxa_ancombc2_basic$res

# 4.3 Process results
bv_res <- res_raw %>%
  dplyr::select(taxon, lfc_GroupBV, q_GroupBV) %>%
  dplyr::filter(abs(lfc_GroupBV) > 1)

bv_res_total <- res_raw_total %>%
  dplyr::select(taxon, lfc_GroupBV, q_GroupBV) %>%
  dplyr::filter(abs(lfc_GroupBV) > 1)

# 4.4 Volcano plot functions
create_volcano_plot <- function(tt_adj, title_suffix = "") {
  top_taxa <- tt_adj %>%
    dplyr::filter(q_GroupBV < 0.05 & abs(lfc_GroupBV) >= 1) %>%
    dplyr::arrange(q_GroupBV) %>%
    head(50)

  tt_adj$Significance <- dplyr::case_when(
    tt_adj$q_GroupBV < 0.001 ~ "p < 0.001",
    tt_adj$q_GroupBV < 0.01 ~ "p < 0.01",
    tt_adj$q_GroupBV < 0.05 ~ "p < 0.05",
    TRUE ~ "Not significant"
  )

  volcano_plot <- ggplot2::ggplot(tt_adj, ggplot2::aes(x = lfc_GroupBV, y = -log10(q_GroupBV))) +
    ggplot2::geom_point(ggplot2::aes(color = Significance), alpha = 0.6, size = 2) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    ggrepel::geom_text_repel(
      data = top_taxa,
      ggplot2::aes(label = taxon),
      max.overlaps = 10,
      box.padding = 0.3,
      point.padding = 0.3,
      size = 3
    ) +
    ggplot2::scale_color_manual(
      values = c("p < 0.001" = "#d73027", "p < 0.01" = "#fc8d59",
                 "p < 0.05" = "#fee08b", "Not significant" = "#999999")
    ) +
    ggplot2::labs(
      x = expression(bold(paste("log"["2"], "(", italic("fold-change"), ")"))),
      y = expression(bold(paste("-log"["10"], "(", italic(p), ")"))),
      color = "Significance"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

  return(volcano_plot)
}

create_volcano_plot_total <- function(tt_adj, title_suffix = "") {
  top_taxa <- tt_adj %>%
    dplyr::filter(q_GroupBV < 0.05 & abs(lfc_GroupBV) >= 1) %>%
    dplyr::arrange(q_GroupBV) %>%
    head(20)

  tt_adj$Significance <- dplyr::case_when(
    tt_adj$q_GroupBV < 0.001 ~ "p < 0.001",
    tt_adj$q_GroupBV < 0.01 ~ "p < 0.01",
    tt_adj$q_GroupBV < 0.05 ~ "p < 0.05",
    TRUE ~ "Not significant"
  )

  volcano_plot <- ggplot2::ggplot(tt_adj, ggplot2::aes(x = lfc_GroupBV, y = -log10(q_GroupBV))) +
    ggplot2::geom_point(ggplot2::aes(color = Significance), alpha = 0.6, size = 2) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    ggrepel::geom_text_repel(
      data = top_taxa,
      ggplot2::aes(label = taxon),
      max.overlaps = 10, box.padding = 0.3, point.padding = 0.3, size = 3
    ) +
    ggplot2::scale_color_manual(
      values = c("p < 0.001" = "#d73027", "p < 0.01" = "#fc8d59",
                 "p < 0.05" = "#fee08b", "Not significant" = "#999999")
    ) +
    ggplot2::labs(
      x = expression(bold(paste("log"["2"], "(", italic("fold-change"), ")"))),
      y = expression(bold(paste("-log"["10"], "(", italic(p), ")"))),
      color = "Significance"
    ) +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      legend.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 16, face = "bold"),
      legend.position = "bottom"
    )

  return(volcano_plot)
}

volcano_bv_pos_taxa_adjusted <- create_volcano_plot(bv_res, "")
ggplot2::ggsave("figures/volcano_da_bv_tax_BV.png",
                volcano_bv_pos_taxa_adjusted, width = 8, height = 6)

volcano_bv_pos_taxa_total <- create_volcano_plot_total(bv_res_total, "")
ggplot2::ggsave("figures/volcano_da_bv_tax_BV_total.png",
                volcano_bv_pos_taxa_total, width = 8, height = 6, dpi = 300)

# 4.5 Export results
res_raw %>%
  dplyr::select(taxon, lfc_GroupBV, q_GroupBV) %>%
  openxlsx::write.xlsx(file = "supplementary_tables/Table_S3.xlsx")

# --- 5. FUNCTIONAL PATHWAY ANALYSIS -----------------------------------------

# 5.1 Load pathway data
load("objects/merged_path_abundance.rda")
load("objects/merged_pathway_coverage.rda")
load("objects/merged_gene_family.rda")

meta_df <- as.data.frame(phyloseq::sample_data(physeq_bact_filtered))

classified_pathways <- merged_path_abundance %>%
  dplyr::filter(!stringr::str_detect(Pathway, "UNMAPPED"),
                !stringr::str_detect(Pathway, "UNINTEGRATED"))

classified_pathways.filtered <- classified_pathways %>%
  dplyr::select(Pathway, dplyr::any_of(rownames(meta_df)))

# 5.2 Create pathway phyloseq object
path_mat <- classified_pathways.filtered %>%
  tibble::column_to_rownames("Pathway") %>%
  as.matrix()

taxonomy <- data.frame(
  Pathway = classified_pathways.filtered$Pathway,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    is_stratified = grepl("\\|", Pathway),
    pathway_id = ifelse(grepl(":", Pathway), sub(":.*$", "", Pathway), Pathway),
    taxonomy = dplyr::case_when(
      is_stratified & stringr::str_detect(Pathway, "s__") ~ "species",
      is_stratified & stringr::str_detect(Pathway, "g__") ~ "genus",
      is_stratified ~ "other_taxonomy",
      TRUE ~ "unstratified"
    ),
    genus = dplyr::case_when(
      stringr::str_detect(Pathway, "g__") ~ stringr::str_extract(Pathway, "(?<=g__)[^.]+"),
      TRUE ~ NA_character_
    ),
    species = dplyr::case_when(
      stringr::str_detect(Pathway, "s__") ~ stringr::str_extract(Pathway, "(?<=s__)[^|]+"),
      TRUE ~ NA_character_
    )
  ) %>%
  tibble::column_to_rownames("Pathway")

pathway_phyloseq <- phyloseq::phyloseq(
  phyloseq::otu_table(path_mat, taxa_are_rows = TRUE),
  phyloseq::sample_data(meta_df),
  phyloseq::tax_table(as.matrix(taxonomy))
)

unstrat_phyloseq_path <- pathway_phyloseq %>%
  phyloseq::subset_taxa(is_stratified == FALSE)

physeq_path_filtered <- filter_phyloseq_path(
  ps = unstrat_phyloseq_path,
  rel_abund_thresh = 0.00005,
  min_samples = 10
)

# 5.3 Load pre-computed ANCOM-BC2 pathway results
load("objects/path_ancombc2_basic.rda")
load("objects/path_ancombc2_adjusted.rda")

# 5.4 Process ANCOM-BC2 pathway results
process_ancombc2_results <- function(ancombc_obj, model_name = "basic") {
  res_df <- ancombc_obj$res %>%
    dplyr::mutate(pathway = taxon, model = model_name) %>%
    dplyr::select(
      pathway, lfc_GroupBV, se_GroupBV, W_GroupBV,
      p_GroupBV, q_GroupBV, diff_GroupBV, passed_ss_GroupBV
    ) %>%
    dplyr::rename(
      log2FC = lfc_GroupBV, SE = se_GroupBV, W_stat = W_GroupBV,
      pval = p_GroupBV, padj = q_GroupBV, is_DA = diff_GroupBV,
      passed_sensitivity = passed_ss_GroupBV
    ) %>%
    dplyr::arrange(padj)

  tax_df <- phyloseq::tax_table(physeq_path_filtered) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("pathway")

  res_annotated <- res_df %>%
    dplyr::left_join(tax_df, by = "pathway") %>%
    dplyr::mutate(
      direction = ifelse(log2FC > 0, "Enriched_in_BV", "Enriched_in_Control"),
      significance = dplyr::case_when(
        padj < 0.001 ~ "***",
        padj < 0.01 ~ "**",
        padj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )

  return(res_annotated)
}

results_basic <- process_ancombc2_results(path_ancombc2_basic, "basic")
results_adjusted <- process_ancombc2_results(path_ancombc2_adjusted, "adjusted")

bv_path_res <- path_ancombc2_adjusted$res %>%
  dplyr::mutate(pathway = taxon) %>%
  dplyr::select(pathway, lfc_GroupBV, se_GroupBV, W_GroupBV,
                p_GroupBV, q_GroupBV, diff_GroupBV, passed_ss_GroupBV) %>%
  dplyr::rename(log2FC = lfc_GroupBV, SE = se_GroupBV, W_stat = W_GroupBV,
                pval = p_GroupBV, padj = q_GroupBV, is_DA = diff_GroupBV,
                passed_sensitivity = passed_ss_GroupBV) %>%
  dplyr::arrange(padj)

# 5.5 Pathway volcano plots
create_pathway_volcano <- function(tt_adj, title_suffix = "") {
  top_path <- tt_adj %>%
    dplyr::filter(padj < 0.05 & abs(log2FC) >= 1) %>%
    dplyr::arrange(padj) %>%
    head(10)

  tt_adj$Significance <- dplyr::case_when(
    tt_adj$padj < 0.001 ~ "p < 0.001",
    tt_adj$padj < 0.01 ~ "p < 0.01",
    tt_adj$padj < 0.05 ~ "p < 0.05",
    TRUE ~ "Not significant"
  )

  volcano_plot <- ggplot2::ggplot(tt_adj, ggplot2::aes(x = log2FC, y = -log10(padj))) +
    ggplot2::geom_point(ggplot2::aes(color = Significance), alpha = 0.6, size = 2) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    ggrepel::geom_text_repel(
      data = top_path, ggplot2::aes(label = pathway),
      max.overlaps = 10, box.padding = 0.3, point.padding = 0.3, size = 3
    ) +
    ggplot2::scale_color_manual(
      values = c("p < 0.001" = "#d73027", "p < 0.01" = "#fc8d59",
                 "p < 0.05" = "#fee08b", "Not significant" = "#999999")
    ) +
    ggplot2::labs(
      x = expression(bold(paste("log"["2"], "(", italic("fold-change"), ")"))),
      y = expression(bold(paste("-log"["10"], "(", italic(p), ")"))),
      color = "Significance"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

  return(volcano_plot)
}

volcano_bv_path_adjusted <- create_pathway_volcano(bv_path_res)
ggplot2::ggsave("figures/volcano_da_path_BV_full.png",
                volcano_bv_path_adjusted, width = 8, height = 6)

# 5.6 Export pathway results
save(bv_path_res, file = "objects/bv_path_res.rda")
openxlsx::write.xlsx(bv_path_res, file = "supplementary_tables/Table_S4.xlsx")

# --- 6. ALPHA DIVERSITY ANALYSIS --------------------------------------------

alpha_div <- phyloseq::estimate_richness(physeq_bact_filtered,
                                         measures = c("Shannon", "Simpson", "Observed"))

alpha_div$Evenness <- alpha_div$Shannon / log(alpha_div$Observed)
alpha_div$Group <- phyloseq::sample_data(physeq_bact_filtered)$Group
alpha_div$Sample_ID <- rownames(alpha_div)

perform_test <- function(metric, group) {
  vals_BV <- metric[group == "BV"]
  vals_Control <- metric[group == "Control"]

  if (length(vals_BV) > 2 && length(vals_Control) > 2) {
    shapiro_BV <- shapiro.test(vals_BV)
    shapiro_Control <- shapiro.test(vals_Control)

    if (shapiro_BV$p.value > 0.05 & shapiro_Control$p.value > 0.05) {
      test_result <- t.test(vals_BV, vals_Control, var.equal = TRUE)
    } else {
      test_result <- wilcox.test(vals_BV, vals_Control)
    }
  } else {
    test_result <- t.test(vals_BV, vals_Control, var.equal = TRUE)
  }

  return(test_result$p.value)
}

p_shannon <- perform_test(alpha_div$Shannon, alpha_div$Group)
p_richness <- perform_test(alpha_div$Observed, alpha_div$Group)
p_evenness <- perform_test(alpha_div$Evenness, alpha_div$Group)

alpha_long <- alpha_div %>%
  dplyr::select(Sample_ID, Group, Shannon, Observed, Evenness) %>%
  tidyr::pivot_longer(cols = c(Shannon, Observed, Evenness),
                      names_to = "Metric",
                      values_to = "Value")

alpha_long$Metric <- factor(alpha_long$Metric,
                            levels = c("Shannon", "Observed", "Evenness"),
                            labels = c("Shannon Diversity", "Richness (Observed)", "Pielou's Evenness"))

p_values <- data.frame(
  Metric = factor(c("Shannon Diversity", "Richness (Observed)", "Pielou's Evenness"),
                  levels = c("Shannon Diversity", "Richness (Observed)", "Pielou's Evenness")),
  p_value = c(p_shannon, p_richness, p_evenness)
)

max_vals <- alpha_long %>%
  dplyr::group_by(Metric) %>%
  dplyr::summarise(max_val = max(Value, na.rm = TRUE)) %>%
  dplyr::left_join(p_values, by = "Metric")

max_vals$y_text <- max_vals$max_val * 1.15
max_vals$y_segment <- max_vals$max_val * 1.10

alpha_plot <- ggplot2::ggplot(alpha_long, ggplot2::aes(x = Group, y = Value, fill = Group)) +
  ggplot2::geom_boxplot() +
  ggplot2::facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::scale_fill_manual(values = c("Control" = "#388ECC", "BV" = "#F68B33")) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 16),
    strip.text = ggplot2::element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.spacing = ggplot2::unit(1, "lines")
  ) +
  ggplot2::geom_text(
    data = max_vals,
    ggplot2::aes(x = 1.5, y = y_text, label = paste("italic(p)*' = '*", sprintf("%.2e", p_value))),
    color = "black", size = 4.5, hjust = 0.5, parse = TRUE, fontface = "bold",
    inherit.aes = FALSE
  ) +
  ggplot2::geom_segment(
    data = max_vals,
    ggplot2::aes(x = 1, xend = 2, y = y_segment, yend = y_segment),
    color = "black", size = 0.4, inherit.aes = FALSE
  )

ggplot2::ggsave("figures/alpha_diversity_panel_plot.png",
                alpha_plot, width = 10, height = 6, dpi = 300)

# --- 7. BETA DIVERSITY ANALYSIS ---------------------------------------------

# 7.1 Bray-Curtis distance
bray_dist <- phyloseq::distance(physeq_bact_filtered, method = "bray")

# 7.2 NMDS ordination
nmds_ord <- phyloseq::ordinate(physeq_bact_filtered, method = "NMDS", distance = bray_dist)

p_beta <- phyloseq::plot_ordination(physeq_bact_filtered, nmds_ord,
                                    type = "samples", color = "Group") +
  ggplot2::geom_point(size = 3) +
  ggplot2::stat_ellipse(level = 0.95, linetype = 2) +
  ggplot2::scale_color_manual(values = c("Control" = "#388ECC", "BV" = "#F68B33")) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = 14),
    axis.title = ggplot2::element_text(size = 16, face = "bold"),
    legend.text = ggplot2::element_text(size = 14),
    legend.title = ggplot2::element_text(size = 16, face = "bold"),
    plot.tag = ggplot2::element_text(size = 18, face = "bold")
  )

# 7.3 PCoA
pcoa_ord <- phyloseq::ordinate(physeq_bact_filtered, method = "PCoA", distance = bray_dist)
var_explained <- pcoa_ord$values$Eigenvalues[1:2] / sum(pcoa_ord$values$Eigenvalues) * 100

p_pcoa <- phyloseq::plot_ordination(physeq_bact_filtered, pcoa_ord,
                                    type = "samples", color = "Group") +
  ggplot2::geom_point(size = 3) +
  ggplot2::stat_ellipse(level = 0.95, linetype = 2) +
  ggplot2::scale_color_manual(values = c("Control" = "#388ECC", "BV" = "#F68B33")) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = 14),
    axis.title = ggplot2::element_text(size = 16, face = "bold"),
    legend.text = ggplot2::element_text(size = 14),
    legend.title = ggplot2::element_text(size = 16, face = "bold"),
    plot.tag = ggplot2::element_text(size = 18, face = "bold")
  ) +
  ggplot2::labs(
    x = paste0("PCoA 1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PCoA 2 (", round(var_explained[2], 1), "%)")
  )

# 7.4 Combined plot
alpha_plot_noleg <- alpha_plot + ggplot2::theme(legend.position = "none")
p_beta_noleg <- p_beta + ggplot2::theme(legend.position = "none")

top_row <- alpha_plot_noleg
bottom_row <- p_beta_noleg + p_pcoa +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")

combined_plot <- (top_row / bottom_row) +
  patchwork::plot_annotation(tag_levels = 'a') &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 16, face = "bold"))

ggplot2::ggsave("figures/Figure1.png",
                combined_plot, width = 12, height = 12, dpi = 300)

# 7.5 PERMANOVA - univariate tests
meta_df_filtered <- meta_df %>%
  dplyr::select(Group, Vaginal_pH, Antibiotic, BV_History, Age)

results_list <- list()
for (var_name in names(meta_df_filtered)) {
  var_data <- meta_df_filtered[[var_name]]
  valid_idx <- !is.na(var_data)

  if (sum(valid_idx) > 0 && length(unique(var_data[valid_idx])) > 1) {
    temp_dist <- as.dist(as.matrix(bray_dist)[valid_idx, valid_idx])
    temp_df <- data.frame(variable = var_data[valid_idx])
    result <- vegan::adonis2(temp_dist ~ variable, data = temp_df, permutations = 999)
    results_list[[var_name]] <- data.frame(
      Variables = var_name,
      P_value = result$`Pr(>F)`[1],
      R2_value = result$R2[1] * 100
    )
  }
}
results_braycurtis_univariate <- do.call(rbind, results_list)
rownames(results_braycurtis_univariate) <- NULL

# 7.6 PERMANOVA - multivariate test
permanova_res <- vegan::adonis2(
  bray_dist ~ Group + Vaginal_pH + Antibiotic + BV_History + Age,
  data = meta_df,
  permutations = 999,
  by = "terms",
  na.action = na.omit
)

results_braycurtis_multivariate_tibble <- tibble::tibble(
  Term = rownames(permanova_res),
  Df = permanova_res$Df,
  SumOfSqs = permanova_res$SumOfSqs,
  R2 = permanova_res$R2 * 100,
  F_stat = permanova_res$F,
  P_value = permanova_res$`Pr(>F)`
)

results_braycurtis_multivariate_tibble %>%
  dplyr::filter(P_value < 0.05 & !is.na(P_value)) %>%
  openxlsx::write.xlsx(file = "supplementary_tables/Table_S2.xlsx")

# --- 8. METABOLOMICS DIFFERENTIAL ANALYSIS ----------------------------------

# 8.1 Load metabolomics data
# Use the public sample metadata (same as loaded earlier)
load("metadata/sample_metadata.rda")
meta <- meta_data

# NOTE: requires raw metabolomics data not included in repository
polar_negative <- data.table::fread('data/polar_negative_results.tsv')
polar_negative.filtered <- polar_negative %>%
  dplyr::filter(User_ID %in% meta$Sample_ID)

# NOTE: requires raw metabolomics data not included in repository
polar_positive <- data.table::fread('data/polar_positive_results.tsv')
polar_positive.filtered <- polar_positive %>%
  dplyr::filter(User_ID %in% meta$Sample_ID)

# 8.2 Normalization functions
non_feature_cols <- c(
  "Platform", "Polarity", "Use_type", "Injection_order", "Sample_type", "Sample_type1",
  "SMMS_ID", "Replicate", "Include", "Extraction_Batch", "Extraction_date", "Tech_Batch",
  "Chrom_Batch", "Norm_Batch", "Analysis_start_date", "Sample_Comments", "Analysis_Comments",
  "Box_ID", "Box_pos", "User_ID", "Volume_of_Sample_uL", "Volume_of_MeOH_uL",
  "Transfered_to_the_vial_uL", "Sample", "Study_Group", "Clinical_Diagnosis", "Analysis_Group"
)

pqn_normalize <- function(X) {
  ref <- matrixStats::colMedians(as.matrix(X), na.rm = TRUE)
  quot <- X / matrix(rep(ref, each = nrow(X)), nrow = nrow(X))
  fac <- matrixStats::rowMedians(as.matrix(quot), na.rm = TRUE)
  X / fac
}

normalize_mode <- function(df_mode, meta, id_col_mode = "User_ID",
                           id_col_meta = "Sample_ID", name_col_meta = "Sample",
                           batch_cols = c("Extraction_Batch", "Tech_Batch"),
                           do_batch_correction = TRUE) {

  df_mode <- as.data.frame(df_mode, stringsAsFactors = FALSE)
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)

  df_mode <- df_mode %>% dplyr::filter(.data$Sample_type %in% c("Sample", "sample"))
  df_mode <- dplyr::select(df_mode, -dplyr::any_of("Sample"))

  feat_cols0 <- setdiff(colnames(df_mode), intersect(non_feature_cols, colnames(df_mode)))
  all_na <- vapply(df_mode %>% dplyr::select(dplyr::all_of(feat_cols0)),
                   function(x) all(is.na(x)), logical(1))
  feat_cols <- feat_cols0[!all_na]

  xwalk <- df_mode %>%
    dplyr::select(dplyr::all_of(
      c(id_col_mode, "Injection_order", batch_cols[batch_cols %in% colnames(df_mode)])
    )) %>%
    dplyr::rename(Sample_ID = dplyr::all_of(id_col_mode)) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(
      meta %>%
        dplyr::select(dplyr::all_of(c(id_col_meta, name_col_meta))) %>%
        dplyr::rename(Sample_ID = dplyr::all_of(id_col_meta),
                      CanonicalSample = dplyr::all_of(name_col_meta)),
      by = "Sample_ID"
    )

  stopifnot(!anyDuplicated(xwalk$CanonicalSample))
  stopifnot(!anyDuplicated(xwalk$Sample_ID))

  meta_mode <- meta %>%
    dplyr::inner_join(
      xwalk %>% dplyr::select(Sample_ID, CanonicalSample,
                              dplyr::any_of(c("Injection_order", batch_cols))),
      by = "Sample_ID"
    ) %>%
    dplyr::mutate(
      dplyr::across(dplyr::any_of(batch_cols), ~factor(.x)),
      Group = factor(Group, levels = c("Control", "BV"))
    )

  M_build <- df_mode %>%
    dplyr::select(dplyr::all_of(c(id_col_mode, feat_cols))) %>%
    dplyr::rename(Sample_ID = dplyr::all_of(id_col_mode)) %>%
    dplyr::inner_join(
      meta_mode %>% dplyr::select(Sample_ID, CanonicalSample),
      by = "Sample_ID",
      relationship = "many-to-one"
    )

  if (!("CanonicalSample" %in% names(M_build))) {
    cand <- intersect(c("CanonicalSample.x", "CanonicalSample.y", "Sample", "Sample.y", "Sample.x"),
                      names(M_build))
    if (length(cand) == 0) stop("Canonical sample column not found after join.")
    M_build <- dplyr::rename(M_build, CanonicalSample = dplyr::all_of(cand[1]))
  }

  M <- M_build %>%
    dplyr::select(CanonicalSample, dplyr::all_of(feat_cols)) %>%
    tibble::column_to_rownames("CanonicalSample") %>%
    as.matrix()

  keep <- colMeans(!is.na(M)) >= 0.7
  if (any(keep)) M <- M[, keep, drop = FALSE]

  M <- pqn_normalize(M)
  M <- log1p(M)
  M <- t(impute::impute.knn(t(M))$data)

  sds <- matrixStats::colSds(M)
  sds[!is.finite(sds) | sds == 0] <- 1
  M <- sweep(M, 2, sqrt(sds), "/")

  if (do_batch_correction) {
    bcols_present <- batch_cols[batch_cols %in% colnames(meta_mode)]
    if (length(bcols_present) > 0) {
      primary_batch <- if ("Extraction_Batch" %in% bcols_present) "Extraction_Batch" else bcols_present[1]

      design_covs <- stats::model.matrix(
        ~ Group + Vaginal_pH + Antibiotic + BV_History + Age,
        data = meta_mode
      )

      M_corrected <- limma::removeBatchEffect(
        x = t(M),
        batch = meta_mode[[primary_batch]],
        design = design_covs
      )

      M <- t(M_corrected)
      message("Batch correction completed using ", primary_batch)
    }
  }

  meta_mode <- meta_mode %>% dplyr::arrange(factor(CanonicalSample, levels = rownames(M)))

  list(
    X_norm = M,
    features = colnames(M),
    samples = rownames(M),
    meta = meta_mode,
    crosswalk = xwalk %>% dplyr::arrange(factor(CanonicalSample, levels = rownames(M)))
  )
}

res_pos <- normalize_mode(polar_positive.filtered, meta)
res_neg <- normalize_mode(polar_negative.filtered, meta)

# 8.3 Positive mode limma analysis
X_pos <- res_pos$X_norm
meta_pos <- res_pos$meta %>%
  dplyr::mutate(
    Group = factor(Group, levels = c("Control", "BV")),
    Antibiotic = factor(Antibiotic, levels = c("No", "Yes")),
    BV_History = factor(BV_History, levels = c("No", "Yes")),
    Vaginal_pH_s = as.numeric(scale(Vaginal_pH)),
    Age_s = as.numeric(scale(Age))
  ) %>%
  dplyr::filter(!is.na(Group)) %>%
  dplyr::arrange(CanonicalSample)

X_pos <- X_pos[meta_pos$CanonicalSample, , drop = FALSE]

design_adj_pos <- stats::model.matrix(~ Group + Vaginal_pH_s + Antibiotic + BV_History + Age_s, data = meta_pos)

fit_adj_pos <- limma::lmFit(base::t(X_pos), design_adj_pos)
fit_adj_pos <- limma::eBayes(fit_adj_pos, robust = TRUE, trend = TRUE)

tt_adj_pos <- limma::topTable(fit_adj_pos, coef = "GroupBV", number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("Metabolite") %>%
  dplyr::mutate(FDR = p.adjust(P.Value, method = "BH"),
                direction = dplyr::if_else(logFC > 0, "Higher in BV", "Higher in Control")) %>%
  dplyr::arrange(FDR)

openxlsx::write.xlsx(tt_adj_pos, file = "supplementary_tables/Table_S5.xlsx")

# 8.4 Negative mode limma analysis
X_neg <- res_neg$X_norm
meta_neg <- res_neg$meta %>%
  dplyr::mutate(
    Group = factor(Group, levels = c("Control", "BV")),
    Antibiotic = factor(Antibiotic, levels = c("No", "Yes")),
    BV_History = factor(BV_History, levels = c("No", "Yes")),
    Vaginal_pH_s = as.numeric(scale(Vaginal_pH)),
    Age_s = as.numeric(scale(Age))
  ) %>%
  dplyr::filter(!is.na(Group)) %>%
  dplyr::arrange(CanonicalSample)

X_neg <- X_neg[meta_neg$CanonicalSample, , drop = FALSE]

design_adj_neg <- stats::model.matrix(~ Group + Vaginal_pH_s + Antibiotic + BV_History + Age_s, data = meta_neg)

fit_adj_neg <- limma::lmFit(base::t(X_neg), design_adj_neg)
fit_adj_neg <- limma::eBayes(fit_adj_neg, robust = TRUE, trend = TRUE)

tt_adj_neg <- limma::topTable(fit_adj_neg, coef = "GroupBV", number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("Metabolite") %>%
  dplyr::mutate(FDR = p.adjust(P.Value, method = "BH"),
                direction = dplyr::if_else(logFC > 0, "Higher in BV", "Higher in Control")) %>%
  dplyr::arrange(FDR)

openxlsx::write.xlsx(tt_adj_neg, file = "supplementary_tables/Table_S6.xlsx")

# 8.5 Metabolomics volcano plots
create_metabolomics_volcano <- function(tt_adj, title_suffix = "") {
  top_metabolites <- tt_adj %>%
    dplyr::filter(adj.P.Val < 0.05 & abs(logFC) >= -1) %>%
    dplyr::arrange(adj.P.Val) %>%
    head(20)

  tt_adj$Significance <- dplyr::case_when(
    tt_adj$adj.P.Val < 0.001 ~ "p < 0.001",
    tt_adj$adj.P.Val < 0.01 ~ "p < 0.01",
    tt_adj$adj.P.Val < 0.05 ~ "p < 0.05",
    TRUE ~ "Not significant"
  )

  volcano_plot <- ggplot2::ggplot(tt_adj, ggplot2::aes(x = logFC, y = -log10(adj.P.Val))) +
    ggplot2::geom_point(ggplot2::aes(color = Significance), alpha = 0.6, size = 2) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    ggrepel::geom_text_repel(
      data = top_metabolites, ggplot2::aes(label = Metabolite),
      max.overlaps = 10, box.padding = 0.3, point.padding = 0.3, size = 3
    ) +
    ggplot2::scale_color_manual(
      values = c("p < 0.001" = "#d73027", "p < 0.01" = "#fc8d59",
                 "p < 0.05" = "#fee08b", "Not significant" = "#999999")
    ) +
    ggplot2::labs(
      x = expression(bold(paste("log"["2"], "(", italic("fold-change"), ")"))),
      y = expression(bold(paste("-log"["10"], "(", italic(p), ")"))),
      color = "Significance"
    ) +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      legend.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 16, face = "bold"),
      legend.position = "bottom"
    )

  return(volcano_plot)
}

volcano_bv_pos_adj <- create_metabolomics_volcano(tt_adj_pos)
volcano_bv_neg_adj <- create_metabolomics_volcano(tt_adj_neg)

ggplot2::ggsave("figures/volcano_da_bv_pos_BV_adj_v2.png",
                volcano_bv_pos_adj, width = 8, height = 6)
ggplot2::ggsave("figures/volcano_da_bv_neg_BV__adj_v2.png",
                volcano_bv_neg_adj, width = 8, height = 6)

combined_volcano_bv_pos_neg <- (volcano_bv_pos_adj + volcano_bv_neg_adj) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = 'A') &
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave("figures/combined_volcano_bv_pos_neg_v2.png",
                combined_volcano_bv_pos_neg, width = 12, height = 6, dpi = 300)

# --- 9. MULTI-OMICS INTEGRATION PREPARATION ---------------------------------

# 9.1 CLR transformation
clr_transform <- function(x) {
  min_detection <- 1 / ncol(x)
  pseudocount <- min_detection * 0.1
  x_pseudo <- x + pseudocount
  gm <- exp(rowMeans(log(x_pseudo)))
  clr_data <- log(x_pseudo / gm)
  return(clr_data)
}

bact_mat <- as.matrix(as.data.frame(phyloseq::otu_table(physeq_bact_filtered)))
taxonomic_clr <- clr_transform(t(bact_mat))

load("objects/merged_path_abundance.rda")
classified_pathways <- merged_path_abundance %>%
  dplyr::filter(!stringr::str_detect(Pathway, "UNMAPPED"),
                !stringr::str_detect(Pathway, "UNINTEGRATED"))

path_mat <- classified_pathways %>%
  dplyr::select(Pathway, dplyr::any_of(rownames(meta_df))) %>%
  tibble::column_to_rownames("Pathway") %>%
  as.matrix()

pathway_clr <- clr_transform(t(path_mat))

# 9.2 Sample alignment
metabolomics_normalized_pos <- X_pos
metabolomics_normalized_neg <- X_neg
meta_data <- meta

common_samples <- intersect(
  intersect(rownames(metabolomics_normalized_neg), rownames(pathway_clr)),
  rownames(taxonomic_clr)
)

cat("Number of common samples:", length(common_samples), "\n")

metabolomics_aligned_pos <- metabolomics_normalized_pos[common_samples, ]
metabolomics_aligned_neg <- metabolomics_normalized_neg[common_samples, ]
pathways_aligned <- pathway_clr[common_samples, ]
taxonomy_aligned <- taxonomic_clr[common_samples, ]

# 9.3 Variance filtering
preprocess_data <- function(data_matrix, var_threshold = 0.01) {
  feature_vars <- apply(data_matrix, 2, var, na.rm = TRUE)
  keep_features <- feature_vars > var_threshold & !is.na(feature_vars)
  cat("Removed", sum(!keep_features), "low variance features out of", ncol(data_matrix), "\n")
  return(data_matrix[, keep_features, drop = FALSE])
}

X <- list(
  metabolomics_pos = preprocess_data(metabolomics_aligned_pos, var_threshold = 0.01),
  metabolomics_neg = preprocess_data(metabolomics_aligned_neg, var_threshold = 0.01),
  pathways = preprocess_data(pathways_aligned, var_threshold = 0.001),
  taxonomy = preprocess_data(taxonomy_aligned, var_threshold = 0.001)
)

sample_groups <- meta_data$Group
names(sample_groups) <- meta_data$Sample
Y_aligned <- sample_groups[common_samples]

save(X, file = "objects/multiomics_views.rda")
save(Y_aligned, file = "objects/group_labels.rda")
save(meta_data, file = "metadata/sample_metadata.rda")

# --- 10. DIABLO DESIGN MATRIX SPECIFICATION ---------------------------------

design_null <- matrix(0, ncol = length(X), nrow = length(X),
                      dimnames = list(names(X), names(X)))
diag(design_null) <- 0

design_full <- matrix(1, ncol = length(X), nrow = length(X),
                      dimnames = list(names(X), names(X)))
diag(design_full) <- 0

design_biological <- matrix(0, ncol = length(X), nrow = length(X),
                            dimnames = list(names(X), names(X)))
design_biological["metabolomics_pos", "metabolomics_neg"] <- 0.1
design_biological["metabolomics_neg", "metabolomics_pos"] <- 0.1
design_biological["taxonomy", "pathways"] <- 1.0
design_biological["pathways", "taxonomy"] <- 1.0
design_biological["pathways", "metabolomics_pos"] <- 0.8
design_biological["metabolomics_pos", "pathways"] <- 0.8
design_biological["pathways", "metabolomics_neg"] <- 0.8
design_biological["metabolomics_neg", "pathways"] <- 0.8
design_biological["taxonomy", "metabolomics_pos"] <- 0.6
design_biological["metabolomics_pos", "taxonomy"] <- 0.6
design_biological["taxonomy", "metabolomics_neg"] <- 0.6
design_biological["metabolomics_neg", "taxonomy"] <- 0.6
diag(design_biological) <- 0

# Data-driven design (empirical correlations)
cat("\n=== COMPUTING EMPIRICAL CORRELATIONS ===\n")
design_empirical <- matrix(NA, ncol = length(X), nrow = length(X),
                           dimnames = list(names(X), names(X)))

for (i in 1:(length(X) - 1)) {
  for (j in (i + 1):length(X)) {
    cat(sprintf("Computing PLS for %s vs %s...\n", names(X)[i], names(X)[j]))
    pls_result <- mixOmics::pls(X[[i]], X[[j]], ncomp = 1)
    cor_value <- cor(pls_result$variates$X, pls_result$variates$Y)
    design_empirical[i, j] <- design_empirical[j, i] <- abs(cor_value)
  }
}

design_datadriven <- design_empirical
design_datadriven[design_datadriven < 0.6] <- 0
diag(design_datadriven) <- 0

# Hybrid design
design_hybrid <- matrix(0, ncol = length(X), nrow = length(X),
                        dimnames = list(names(X), names(X)))
design_hybrid["metabolomics_pos", "metabolomics_neg"] <- 0.7
design_hybrid["metabolomics_neg", "metabolomics_pos"] <- 0.7
design_hybrid["taxonomy", "pathways"] <- 1.0
design_hybrid["pathways", "taxonomy"] <- 1.0
design_hybrid["pathways", "metabolomics_pos"] <- 0.9
design_hybrid["metabolomics_pos", "pathways"] <- 0.9
design_hybrid["pathways", "metabolomics_neg"] <- 0.9
design_hybrid["metabolomics_neg", "pathways"] <- 0.9
design_hybrid["taxonomy", "metabolomics_pos"] <- 0.85
design_hybrid["metabolomics_pos", "taxonomy"] <- 0.85
design_hybrid["taxonomy", "metabolomics_neg"] <- 0.84
design_hybrid["metabolomics_neg", "taxonomy"] <- 0.84
diag(design_hybrid) <- 0

design_list <- list(
  null = design_null, full = design_full, biological = design_biological,
  datadriven = design_datadriven, hybrid = design_hybrid
)

# --- 11. DIABLO MODEL FITTING AND TUNING ------------------------------------

load(file = "objects/tune_ncomp_full_design_BV_v3.rda")
tune_full <- tune_result
load(file = "objects/tune_ncomp_biological_design_BV_v3.rda")
tune_biological <- tune_result
load(file = "objects/tune_ncomp_hybrid_design_BV_v3.rda")
tune_hybrid <- tune_result
load(file = "objects/tune_ncomp_datadriven_design_BV_v3.rda")
tune_datadriven <- tune_result
load(file = "objects/tune_ncomp_null_design_BV_v3.rda")
tune_null <- tune_result

# --- 12. DIABLO PERFORMANCE EVALUATION --------------------------------------

comparison_df <- data.frame(
  Design = c("null", "full", "biological", "datadriven", "hybrid"),
  Optimal_ncomp = c(
    tune_null$choice.ncomp$ncomp, tune_full$choice.ncomp$ncomp,
    tune_biological$choice.ncomp$ncomp, tune_datadriven$choice.ncomp$ncomp,
    tune_hybrid$choice.ncomp$ncomp
  ),
  BER_comp1 = c(tune_null$error.rate[1], tune_full$error.rate[1],
                tune_biological$error.rate[1], tune_datadriven$error.rate[1],
                tune_hybrid$error.rate[1]),
  BER_comp2 = c(tune_null$error.rate[2], tune_full$error.rate[2],
                tune_biological$error.rate[2], tune_datadriven$error.rate[2],
                tune_hybrid$error.rate[2]),
  BER_comp3 = c(tune_null$error.rate[3], tune_full$error.rate[3],
                tune_biological$error.rate[3], tune_datadriven$error.rate[3],
                tune_hybrid$error.rate[3])
)

comparison_df <- comparison_df[order(comparison_df$BER_comp2), ]
print(comparison_df, row.names = FALSE)

# --- 13. DIABLO VISUALIZATION -----------------------------------------------

# 13.1 Fit final model with optimal parameters
names(X) <- c("Metab_pos", "Metab_neg", "Path", "Taxa")

final_ncomp <- 3
final_keepX <- list(
  Metab_pos = c(30, 5, 30),
  Metab_neg = c(5, 20, 30),
  Path      = c(10, 5, 20),
  Taxa      = c(5, 5, 20)
)

final_diablo_model <- mixOmics::block.splsda(
  X = X, Y = Y_aligned,
  ncomp = final_ncomp,
  keepX = final_keepX,
  design = design_datadriven
)

# 13.2 Performance assessment
perf_final <- mixOmics::perf(
  final_diablo_model,
  validation = "Mfold",
  folds = 5,
  nrepeat = 50,
  cpus = 4
)

cat("\n=== FINAL MODEL PERFORMANCE ===\n")
print(perf_final$error.rate)

# 13.3 Feature extraction
selected_features_comp1 <- mixOmics::selectVar(final_diablo_model, comp = 1)
selected_features_comp2 <- mixOmics::selectVar(final_diablo_model, comp = 2)
selected_features_comp3 <- mixOmics::selectVar(final_diablo_model, comp = 3)

extract_block_features <- function(selected_obj, block_name, comp_num) {
  if (!is.null(selected_obj[[block_name]])) {
    selected_obj[[block_name]]$value %>%
      tibble::rownames_to_column("Feature") %>%
      dplyr::mutate(Component = comp_num, Block = block_name,
                    Abs_Loading = abs(value.var)) %>%
      dplyr::arrange(dplyr::desc(Abs_Loading))
  } else { NULL }
}

all_selected_features <- dplyr::bind_rows(
  extract_block_features(selected_features_comp1, "Metab_pos", 1),
  extract_block_features(selected_features_comp1, "Metab_neg", 1),
  extract_block_features(selected_features_comp1, "Path", 1),
  extract_block_features(selected_features_comp1, "Taxa", 1),
  extract_block_features(selected_features_comp2, "Metab_pos", 2),
  extract_block_features(selected_features_comp2, "Metab_neg", 2),
  extract_block_features(selected_features_comp2, "Path", 2),
  extract_block_features(selected_features_comp2, "Taxa", 2),
  extract_block_features(selected_features_comp3, "Metab_pos", 3),
  extract_block_features(selected_features_comp3, "Metab_neg", 3),
  extract_block_features(selected_features_comp3, "Path", 3),
  extract_block_features(selected_features_comp3, "Taxa", 3)
)

openxlsx::write.xlsx(all_selected_features, file = "supplementary_tables/Table_S8.xlsx")

# Variance explained
var_explained_diablo <- final_diablo_model$prop_expl_var
var_df <- purrr::map2_dfr(var_explained_diablo, names(var_explained_diablo),
                           ~ data.frame(Block = .y,
                                        Component = paste0("Comp", 1:length(.x)),
                                        Variance_Explained = .x))
var_df %>%
  dplyr::filter(!Block == "Y") %>%
  dplyr::mutate(Variance_Explained = round(Variance_Explained * 100, 1)) %>%
  openxlsx::write.xlsx(file = "supplementary_tables/Table_S7.xlsx")

# 13.4 Sample plots
png("figures/integration_plots/DIABLO_sample_plot_comp1-2.png",
    width = 1200, height = 1200, res = 150)
mixOmics::plotIndiv(final_diablo_model, comp = c(1, 2), group = Y_aligned,
                    ind.names = FALSE, legend = TRUE,
                    title = "DIABLO Multi-Omics Integration: BV vs Control",
                    style = "ggplot2", size.title = 14)
dev.off()

# 13.5 AUROC analysis
auroc_results <- mixOmics::auroc(final_diablo_model)

auroc_data <- data.frame()
blocks <- c("Metab_pos", "Metab_neg", "Path", "Taxa")
block_labels <- c("Metabolomics\n(Positive)", "Metabolomics\n(Negative)",
                  "Pathways", "Taxonomy")

for (i in 1:length(blocks)) {
  block <- blocks[i]
  block_label <- block_labels[i]
  for (comp in 1:3) {
    comp_name <- paste0("comp", comp)
    if (!is.null(auroc_results[[block]][[comp_name]])) {
      auc <- auroc_results[[block]][[comp_name]]["Control vs BV", "AUC"]
      pval <- auroc_results[[block]][[comp_name]]["Control vs BV", "p-value"]
      auroc_data <- rbind(auroc_data, data.frame(
        Block = block_label, Block_short = block, Component = comp,
        AUC = auc, P_value = pval,
        Sig = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**",
                     ifelse(pval < 0.05, "*", "ns")))
      ))
    }
  }
}
auroc_data$Performance <- cut(auroc_data$AUC, breaks = c(0, 0.7, 0.8, 0.9, 1.0),
                              labels = c("Poor", "Fair", "Good", "Excellent"))

# AUROC 4-panel figure
block_colors <- c("Metabolomics\n(Positive)" = "#2E86AB", "Metabolomics\n(Negative)" = "#A23B72",
                  "Pathways" = "#F18F01", "Taxonomy" = "#06A77D")

plot_auroc_A <- ggplot2::ggplot(auroc_data, ggplot2::aes(x = factor(Component), y = AUC, fill = Block)) +
  ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(0.9), color = "black", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", AUC)),
                     position = ggplot2::position_dodge(0.9), vjust = -0.5, size = 3, fontface = "bold") +
  ggplot2::geom_hline(yintercept = 0.9, linetype = "dashed", color = "red", linewidth = 0.8) +
  ggplot2::scale_fill_manual(values = block_colors) +
  ggplot2::scale_y_continuous(limits = c(0.85, 1.0), breaks = seq(0.85, 1.0, 0.05)) +
  ggplot2::labs(title = "A. Block-Specific AUROC Performance", x = "Component",
                y = "Area Under ROC Curve (AUC)", fill = "Data Block") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14), legend.position = "right")

plot_auroc_B <- ggplot2::ggplot(auroc_data, ggplot2::aes(x = factor(Component), y = Block, fill = AUC)) +
  ggplot2::geom_tile(color = "white", linewidth = 1) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f\n%s", AUC, Sig)),
                     size = 3.5, fontface = "bold", color = "white") +
  viridis::scale_fill_viridis(option = "plasma", begin = 0.3, end = 0.9, name = "AUC") +
  ggplot2::labs(title = "B. AUC Heatmap Across Components", x = "Component", y = "Data Block") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14), panel.grid = ggplot2::element_blank())

pdf("figures/integration_plots/DIABLO_AUROC_comprehensive.pdf", width = 16, height = 12)
gridExtra::grid.arrange(plot_auroc_A, plot_auroc_B, ncol = 2,
                        top = grid::textGrob("DIABLO Model: Block-Specific AUROC Performance",
                                             gp = grid::gpar(fontsize = 18, fontface = "bold")))
dev.off()

# 13.6 Circos plots
pdf("figures/integration_plots/circos_plot_diablo_v4.pdf", width = 10, height = 10)
mixOmics::circosPlot(final_diablo_model, cutoff = 0.8, line = TRUE, size.variables = 0.6)
dev.off()

png("figures/integration_plots/circos_plot_diablo_v3.png", width = 800, height = 800, res = 100)
mixOmics::circosPlot(final_diablo_model, cutoff = 0.8, line = TRUE, size.variables = 0.6)
dev.off()

# CIM heatmap
mixOmics::cimDiablo(final_diablo_model, cutoff = 0.8)

# Block-specific circos plots
for (cutoff_val in c(0.5, 0.55)) {
  block_pairs <- list(c("Taxa", "Metab_pos"), c("Taxa", "Metab_neg"),
                      c("Taxa", "Path"), c("Path", "Metab_pos"), c("Path", "Metab_neg"))
  for (bp in block_pairs) {
    fname <- paste0("figures/integration_plots/circos_plot_cutoff_",
                    cutoff_val, "_", tolower(bp[1]), "_", tolower(bp[2]), "_comp1-2-3.png")
    png(fname, width = 1200, height = 1200, res = 150)
    mixOmics::circosPlot(final_diablo_model, cutoff = cutoff_val,
                         blocks = bp, line = TRUE, size.variables = 0.6,
                         legend.title = "Abundance", showIntraLinks = FALSE)
    dev.off()
  }
}

# 13.7 Feature importance extraction
extract_top_features_direct <- function(model, n_features = 15, components = c(1, 2)) {
  top_features <- list()
  for (block in names(model$X)) {
    loadings_matrix <- model$loadings[[block]]
    feature_importance <- apply(abs(loadings_matrix[, components, drop = FALSE]), 1, sum)
    top_idx <- order(feature_importance, decreasing = TRUE)[1:min(n_features, length(feature_importance))]
    top_features[[block]] <- rownames(loadings_matrix)[top_idx]
    attr(top_features[[block]], "importance") <- feature_importance[top_idx]
  }
  return(top_features)
}

cat("\n=== TOP FEATURES ===\n")
print(extract_top_features_direct(final_diablo_model, n_features = 20, components = c(1, 2)))

# --- 14. GARDNERELLA PANGENOME ANALYSIS --------------------------------------
# Moved to scripts/02_genome_metabolome.R
# Outputs: supplementary_tables/Table_S11.xlsx, Figure 7
# -----------------------------------------------------------------------------

# --- 15. SESSION INFO -------------------------------------------------------

cat("\n\n=== SESSION INFO ===\n")
sessionInfo()

cat("\n\n=== ANALYSIS COMPLETE ===\n")
cat("All figures saved to: figures/\n")
cat("All tables saved to: supplementary_tables/\n")
cat("All objects saved to: objects/ and metadata/\n")
