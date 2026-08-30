#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# QA LAYER 2b - VALUE-LEVEL LEAK DETECTION (author-side)
#
# Name-based scanning (scan_rda.R) cannot catch a withheld variable whose column
# has been renamed, or one hidden in an attribute. This layer takes the ACTUAL
# withheld vectors from the private files and searches every publishable object
# for them numerically: exact match on the sorted values, or a perfect linear
# rescale (which catches z-scored or unit-converted copies).
#
# It requires the private inputs, so it runs only on an author's machine and is
# skipped elsewhere (including in CI, which never has them).
#
# Usage:  Rscript scripts/qa/scan_values.R
# Exit 0 = clean or skipped, 1 = leak.
# ---------------------------------------------------------------------------

PRIV_META   <- "metadata/sample_metadata.rda"
PRIV_GENOME <- "objects/genome_metabolome_results.rda"

if (!file.exists(PRIV_META) && !file.exists(PRIV_GENOME)) {
  cat("QA L2b - value scan: SKIPPED (private inputs absent; nothing to compare against)\n")
  quit(status = 0)
}

targets <- list()
add <- function(name, v) {
  v <- suppressWarnings(as.numeric(v)); v <- v[!is.na(v)]
  if (length(v) >= 20) targets[[name]] <<- sort(v)
}
if (file.exists(PRIV_META)) {
  e <- new.env(); load(PRIV_META, envir = e); m <- as.data.frame(get("meta_data", envir = e))
  for (cn in setdiff(colnames(m), c("Sample", "Group"))) add(paste0("meta$", cn), m[[cn]])
}
if (file.exists(PRIV_GENOME)) {
  g <- new.env(); load(PRIV_GENOME, envir = g)
  if ("covariates" %in% ls(g)) {
    cv <- as.data.frame(get("covariates", envir = g))
    for (cn in setdiff(colnames(cv), c("Sample", "Group", "MAG_Name", "Species", "Base_Sample_ID")))
      add(paste0("covariates$", cn), cv[[cn]])
  }
}

pub <- suppressWarnings(system2("git", c("ls-files", "--cached", "--others", "--exclude-standard"),
                                stdout = TRUE, stderr = FALSE))
pub <- grep("[.](rda|RData|rds)$", pub, value = TRUE, ignore.case = TRUE)
pub <- setdiff(pub, c(PRIV_META, PRIV_GENOME))

vecs <- list()
collect <- function(o, path, d = 0L) {
  if (d > 12L) return(invisible(NULL))
  if (isS4(o)) {
    for (sn in methods::slotNames(o))
      try(collect(methods::slot(o, sn), paste0(path, "@", sn), d + 1L), silent = TRUE)
    return(invisible(NULL))
  }
  if (is.environment(o)) {
    # A formula, lm or ggplot saved from a scope holding the clinical table
    # carries it in .Environment; without this branch that is invisible.
    # globalenv/baseenv/namespaces are session state and must not be walked.
    if (identical(o, globalenv()) || identical(o, baseenv()) ||
        identical(o, emptyenv()) || environmentName(o) != "" || isNamespace(o))
      return(invisible(NULL))
    for (n in ls(o, all.names = TRUE))
      try(collect(get(n, envir = o), paste0(path, "<env>$", n), d + 1L), silent = TRUE)
    return(invisible(NULL))
  }
  if (is.matrix(o) && is.numeric(o)) {
    # A matrix must be split: flattening it hides a renamed clinical column.
    cn <- colnames(o); rn <- rownames(o)
    if (ncol(o) >= 1 && nrow(o) >= 20)
      for (j in seq_len(ncol(o)))
        vecs[[sprintf("%s[,%s]", path, if (!is.null(cn)) cn[j] else j)]] <<- as.numeric(o[, j])
    if (nrow(o) >= 1 && ncol(o) >= 20)
      for (i in seq_len(nrow(o)))
        vecs[[sprintf("%s[%s,]", path, if (!is.null(rn)) rn[i] else i)]] <<- as.numeric(o[i, ])
  }
  if (is.numeric(o) && !is.matrix(o) && length(o) >= 20) vecs[[path]] <<- as.numeric(o)
  if (is.data.frame(o)) for (cn in colnames(o)) {
    v <- o[[cn]]; if (is.numeric(v) && length(v) >= 20) vecs[[paste0(path, "$", cn)]] <<- as.numeric(v)
  }
  # Attributes can carry data too; a bare vector with an `aux` attribute was a
  # documented blind spot before this branch existed.
  for (an in names(attributes(o)))
    if (!an %in% c("names", "class", "row.names", "dim", "dimnames", "levels"))
      try(collect(attr(o, an), paste0(path, "@attr:", an), d + 1L), silent = TRUE)
  if (is.list(o)) {
    nm <- names(o); if (is.null(nm)) nm <- paste0("[[", seq_along(o), "]]")
    for (i in seq_along(o)) try(collect(o[[i]], paste0(path, "$", nm[i]), d + 1L), silent = TRUE)
  }
  invisible(NULL)
}
for (f in pub) {
  en <- new.env()
  ok <- try(if (grepl("[.]rds$", f, ignore.case = TRUE)) assign("x", readRDS(f), envir = en)
            else load(f, envir = en), silent = TRUE)
  if (inherits(ok, "try-error")) next
  for (n in ls(en)) try(collect(get(n, envir = en), sprintf("%s :: %s", f, n)), silent = TRUE)
}

cat(sprintf("QA L2b - value scan: %d withheld vector(s) vs %d numeric vector(s) in %d publishable object(s)\n",
            length(targets), length(vecs), length(pub)))

hits <- character(0)
for (tn in names(targets)) {
  tv <- targets[[tn]]
  for (vn in names(vecs)) {
    v <- vecs[[vn]]; v <- v[!is.na(v)]
    if (length(v) != length(tv)) next
    sv <- sort(v)
    exact <- isTRUE(all.equal(sv, tv, tolerance = 1e-8))
    resc  <- FALSE
    if (!exact && stats::sd(sv) > 0 && stats::sd(tv) > 0) {
      r <- suppressWarnings(abs(stats::cor(sv, tv)))
      resc <- !is.na(r) && r > 0.9999
    }
    if (exact || resc)
      hits <- c(hits, sprintf("  LEAK  %s  %s  %s (n=%d)", vn,
                              if (exact) "==" else "~=", tn, length(tv)))
  }
}
if (length(hits)) {
  cat("FAIL - withheld values found in publishable objects:\n")
  cat(paste(unique(hits), collapse = "\n"), "\n")
  quit(status = 1)
}
cat("PASS - no withheld values present, exact or rescaled.\n")
quit(status = 0)
