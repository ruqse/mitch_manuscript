# ---------------------------------------------------------------------------
# _metadata_access.R  -  shared metadata / cached-object loader
#
# This repository ships DE-IDENTIFIED objects. Individual-level clinical
# variables (Age, Vaginal_pH, Antibiotic, BV_History) were not published with
# the article and are available only under controlled access -- see
# "Controlled access" in README.md.
#
# Each loader prefers the private participant-level file when it is present
# (i.e. on an author's machine) and otherwise falls back to the public
# derivative that is committed here. Analyses that require the withheld
# variables are gated on HAS_CLINICAL rather than failing.
#
# Sourced by scripts/01_MOFA2_integration.R and scripts/02_genome_metabolome.R
# ---------------------------------------------------------------------------

CLINICAL_VARS <- c("Age", "Vaginal_pH", "Antibiotic", "BV_History")

.pick <- function(private, public, what) {
  if (file.exists(private)) return(private)
  if (file.exists(public))  return(public)
  stop(sprintf("Neither '%s' nor '%s' found - cannot load %s.", private, public, what),
       call. = FALSE)
}

load_sample_metadata <- function() {
  f <- .pick("metadata/sample_metadata.rda",
             "metadata/sample_metadata_public.rda", "sample metadata")
  e <- new.env(); load(f, envir = e)
  md <- as.data.frame(get("meta_data", envir = e), stringsAsFactors = FALSE)
  message(sprintf("[metadata] %s  (%d x %d: %s)",
                  f, nrow(md), ncol(md), paste(colnames(md), collapse = ", ")))
  md
}

load_genome_metabolome_cache <- function(envir = parent.frame()) {
  f <- .pick("objects/genome_metabolome_results.rda",
             "objects/genome_metabolome_results_public.rda",
             "genome-metabolome results")
  load(f, envir = envir)
  message(sprintf("[cache] %s", f))
  invisible(f)
}

has_clinical <- function(md) all(CLINICAL_VARS %in% colnames(md))

skip_clinical <- function(what) {
  message(sprintf(
    paste0("\n[SKIPPED] %s\n",
           "  Requires participant-level clinical variables (%s), which are not\n",
           "  distributed in this repository. They are available under controlled\n",
           "  access from the corresponding author; see README.md, 'Controlled access'.\n"),
    what, paste(CLINICAL_VARS, collapse = ", ")))
  invisible(NULL)
}
