#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# make_public_objects.R
#
# Derives the PUBLIC, de-identified R objects that this repository ships from
# the private, participant-level objects held only on the authors' machines.
#
#   INPUT  (never committed; see .gitignore)
#     metadata/sample_metadata.rda            meta_data    111 x 7
#     objects/genome_metabolome_results.rda   15 objects, incl. `covariates` 33 x 13
#
#   OUTPUT (committed)
#     metadata/sample_metadata_public.rda
#     objects/genome_metabolome_results_public.rda
#
# WHAT IS REMOVED AND WHY
# -----------------------
# The study was approved by the Swedish Ethical Review Authority
# (Etikprovningsmyndigheten), Lund Division 2 Medicine, approval 2021-05794-01.
# Individual-level clinical variables were NOT published with the article --
# they appear there only as model covariates and aggregate results. Publishing
# them here, keyed to sample identifiers that also key the ENA submission
# (PRJEB108308), would be a new disclosure of special-category health data.
# They are therefore withheld and available only under controlled access.
#
# Removed variables: Sample_ID, Base_Sample_ID, Vaginal_pH, Age, Antibiotic,
#                    BV_History, pH_z, Age_z, Group_BV
#
# Retained: Sample, Group, MAG_Name, Species.
#
# On `Group` (BV/Control), stated precisely: the article's Supplementary Data 3 /
# Table S4 already publishes Sample -> BV_status, but only for the 53 samples
# carrying Gardnerella MAGs (35 in-cohort MAGs, from 33 distinct samples) -- NOT
# for all 111. So the
# 111-row Sample -> Group table shipped here is broader than that precedent.
# It is retained because Group is the study's primary grouping variable, is
# recomputable from the taxonomy matrix shipped in objects/multiomics_views.rda
# and from the public ENA metagenomes, and attaches only to a pseudonymous
# sample identifier. It carries no clinical measurement.
#
# NOTE ON dbrda_cond
# ------------------
# `dbrda_cond` is a vegan::capscale object fitted with
#   Condition(pH_z + Age_z + Antibiotic + BV_History).
# Its $pCCA$QR component stores the QR decomposition of that conditioning
# matrix (33 x 4, columns pH_z/Age_z/Antibiotic/BV_History). A QR
# decomposition permits reconstruction of the original values, so shipping the
# fitted object would leak the very variables withheld above. The object is
# therefore DROPPED and replaced by the single scalar the analysis actually
# uses (its adjusted R-squared), alongside the permutation test `anova_cond`.
#
# `dbrda_uncond` is unconditioned ($pCCA is NULL) and is safe.
# `varpart_result` stores only variance fractions and the call text - no
# per-sample values - and is safe.
#
# Usage:  Rscript scripts/make_public_objects.R
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages(library(vegan)))

IN_META   <- "metadata/sample_metadata.rda"
IN_GENOME <- "objects/genome_metabolome_results.rda"
OUT_META   <- "metadata/sample_metadata_public.rda"
OUT_GENOME <- "objects/genome_metabolome_results_public.rda"

for (f in c(IN_META, IN_GENOME)) {
  if (!file.exists(f)) {
    stop(sprintf(
      "Missing private input '%s'.\nThis script is for the authors only; the public repository ships its outputs, not its inputs.",
      f), call. = FALSE)
  }
}

PUBLIC_META_COLS       <- c("Sample", "Group")
PUBLIC_COVARIATE_COLS  <- c("Sample", "MAG_Name", "Species", "Group")

# --- 1. sample metadata -----------------------------------------------------
e <- new.env(); load(IN_META, envir = e)
meta_data <- get("meta_data", envir = e)[, PUBLIC_META_COLS]
meta_data <- as.data.frame(meta_data, stringsAsFactors = FALSE)
cat(sprintf("sample_metadata: %d x %d -> %d x %d (%s)\n",
            nrow(get("meta_data", envir = e)), ncol(get("meta_data", envir = e)),
            nrow(meta_data), ncol(meta_data),
            paste(colnames(meta_data), collapse = ", ")))
save(meta_data, file = OUT_META)

# --- 2. genome-metabolome results ------------------------------------------
g <- new.env(); load(IN_GENOME, envir = g)

covariates <- as.data.frame(get("covariates", envir = g))[, PUBLIC_COVARIATE_COLS]
cat(sprintf("covariates:      33 x 13 -> %d x %d (%s)\n",
            nrow(covariates), ncol(covariates),
            paste(colnames(covariates), collapse = ", ")))

# Preserve the scalar that dbrda_cond was used for, then discard the object.
dbrda_cond_adjR2 <- RsquareAdj(get("dbrda_cond", envir = g))$adj.r.squared
cat(sprintf("dbrda_cond:      DROPPED (leaks conditioning matrix via $pCCA$QR); adj R2 %.6f retained\n",
            dbrda_cond_adjR2))

keep <- setdiff(ls(g), c("covariates", "dbrda_cond"))
for (n in keep) assign(n, get(n, envir = g))

save(list = c(keep, "covariates", "dbrda_cond_adjR2"), file = OUT_GENOME)
cat(sprintf("wrote %s (%d objects)\n", OUT_GENOME, length(keep) + 2))
cat("\nDone. Run scripts/qa/scan_rda.R to verify no clinical columns remain.\n")
