#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# QA LAYER 2 - deep binary object scan
#
# Recursively walks EVERY nested object inside every .rda / .rds / .xlsx that
# git would publish, and fails if a participant-level clinical variable name
# appears as a column, name, or dimname at ANY depth.
#
# This layer exists because .gitignore guards paths, not contents. Two real
# leaks in this repository were only visible from inside binary objects:
#   * `covariates` (33 x 13 clinical table) nested in a *results* .rda
#   * `dbrda_cond$pCCA$QR$qr`, a 33 x 4 QR decomposition whose column names
#     are pH_z/Age_z/Antibiotic/BV_History and from which the raw conditioning
#     values can be reconstructed
#
# Usage:  Rscript scripts/qa/scan_rda.R          # scans what git would publish
#         Rscript scripts/qa/scan_rda.R --all    # scans every file on disk
# Exit 0 = clean, 1 = leak found.
# ---------------------------------------------------------------------------

# Limited to variables the published article itself names (SI Table 6:
# "~ Group + vaginal pH + antibiotic + BV history + age") plus generic
# identifiers. Cohort variables the article does not mention are guarded in the
# private authoring repository's copy, so this public file does not enumerate them.
# `group_bv` is deliberately NOT listed: it is a 0/1 coding of Group, and the
# BV/Control label is published (Supplementary Data 3 / Table S4), so it
# discloses nothing further. Listing it also false-positives on the published
# ANCOM-BC2 column names log2FC_GroupBV and q_GroupBV.
BANNED <- c("vaginal_ph", "ph_z", "age", "age_z", "antibiotic",
            "bv_history", "sample_id", "base_sample_id",
            "diagnosis", "personnummer", "dob", "date_of_birth",
            "ethnicity", "parity", "gravidity", "bmi", "hiv", "pregnan")

# `Group`/`BV_status` are intentionally allowed: BV status per sample is already
# public via the article's Supplementary Data 3 / Table S4.

# One alternation, compiled once. Word-anchored so "Average" does not match
# "age", and "_" is optional so "Vaginal pH" matches "vaginal_ph".
BANNED_RE <- paste0("(^|[^a-z0-9])(",
                    paste(gsub("_", "[ _]?", BANNED), collapse = "|"),
                    ")([^a-z0-9]|$)")

# Environments that must never be traversed: doing so walks the whole session.
.is_global_env <- function(e) {
  identical(e, globalenv()) || identical(e, baseenv()) || identical(e, emptyenv()) ||
    environmentName(e) != "" || isNamespace(e)
}

args    <- commandArgs(trailingOnly = TRUE)
scan_all <- "--all" %in% args

files <- if (scan_all) {
  list.files(".", pattern = "[.](rda|RData|rds|xlsx)$", recursive = TRUE,
             ignore.case = TRUE, full.names = FALSE)
} else {
  out <- suppressWarnings(system2("git", c("ls-files", "--cached", "--others",
                                           "--exclude-standard"),
                                  stdout = TRUE, stderr = FALSE))
  grep("[.](rda|RData|rds|xlsx)$", out, value = TRUE, ignore.case = TRUE)
}
files <- files[!grepl("^[.]git/", files)]

findings <- list()

check_names <- function(nms, path, kind) {
  if (is.null(nms)) return(invisible(NULL))
  low <- tolower(trimws(nms))
  hit <- nms[low %in% BANNED | grepl(BANNED_RE, low, perl = TRUE)]
  if (length(hit)) {
    findings[[length(findings) + 1L]] <<-
      sprintf("  LEAK  %s  [%s] -> %s", path, kind, paste(unique(hit), collapse = ", "))
  }
}

walk <- function(o, path, depth = 0L) {
  if (depth > 12L) return(invisible(NULL))
  cls <- class(o)[1]

  if (isS4(o)) {
    for (sn in methods::slotNames(o)) {
      try(walk(methods::slot(o, sn), paste0(path, "@", sn), depth + 1L), silent = TRUE)
    }
    return(invisible(NULL))
  }

  if (is.data.frame(o) || is.matrix(o)) {
    check_names(colnames(o), path, paste0(cls, " ", paste(dim(o), collapse = "x")))
    check_names(rownames(o), path, paste0(cls, " rownames"))
  }
  if (!is.null(names(o))) check_names(names(o), path, cls)
  if (!is.null(dim(o)) && !is.null(dimnames(o))) {
    for (d in dimnames(o)) check_names(d, path, paste0(cls, " dimnames"))
  }

  for (an in names(attributes(o)))
    if (!an %in% c("names", "class", "row.names", "dim", "dimnames", "levels"))
      try(walk(attr(o, an), paste0(path, "@attr:", an), depth + 1L), silent = TRUE)

  if (is.environment(o)) {
    # Only local/captured environments are of interest. globalenv, baseenv and
    # namespaces are session state, not repository content.
    if (.is_global_env(o)) return(invisible(NULL))
    for (n in ls(o, all.names = TRUE))
      try(walk(get(n, envir = o), paste0(path, "<env>$", n), depth + 1L), silent = TRUE)
    return(invisible(NULL))
  }

  if (is.list(o)) {
    nms <- names(o)
    if (is.null(nms)) nms <- paste0("[[", seq_along(o), "]]")
    for (i in seq_along(o)) {
      try(walk(o[[i]], paste0(path, "$", nms[i]), depth + 1L), silent = TRUE)
    }
  }
  invisible(NULL)
}

cat("QA L2 - deep object scan\n")
cat(sprintf("scanning %d binary file(s)%s\n\n", length(files),
            if (scan_all) " (all on disk)" else " git would publish"))

for (f in files) {
  if (grepl("[.]xlsx$", f, ignore.case = TRUE)) {
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      for (sh in openxlsx::getSheetNames(f)) {
        d <- try(openxlsx::read.xlsx(f, sheet = sh, rows = 1:3), silent = TRUE)
        if (!inherits(d, "try-error") && !is.null(d)) {
          check_names(colnames(d), sprintf("%s[%s]", f, sh), "xlsx header")
          check_names(as.character(unlist(d[1, , drop = TRUE])),
                      sprintf("%s[%s]", f, sh), "xlsx row1")
        }
      }
    }
    next
  }
  e <- new.env()
  ok <- try(if (grepl("[.]rds$", f, ignore.case = TRUE))
              assign("x", readRDS(f), envir = e) else load(f, envir = e),
            silent = TRUE)
  if (inherits(ok, "try-error")) { cat(sprintf("  SKIP  %s (unreadable)\n", f)); next }
  for (n in ls(e)) try(walk(get(n, envir = e), sprintf("%s :: %s", f, n)), silent = TRUE)
}

if (length(findings)) {
  cat("FAIL - participant-level variables found in publishable objects:\n\n")
  cat(paste(unlist(findings), collapse = "\n"), "\n\n")
  quit(status = 1)
}
cat("PASS - no participant-level clinical variables in any publishable object.\n")
quit(status = 0)
