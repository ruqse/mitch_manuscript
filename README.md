# MiTCH: Microbiome of ITCHing

Multi-omics integration of metagenomic, metabolomic and *Gardnerella* pangenomic
data in bacterial vaginosis (BV). Analysis code and processed data for the MiTCH
study.

## Publication

> Shabana, H.\*, Dube, F.\*, Lahtinen, E., Bergman Røthe, E., Engstrand, L.,
> Warchavchik Hugerth, L., Schuppe-Koistinen, I.<sup>†</sup> &
> Rodriguez-Wallberg, K. A.<sup>†</sup>
> **Multi-omics investigation of metabolomic and microbial features in bacterial
> vaginosis.** *npj Women's Health* **4**, 26 (2026).
> [doi:10.1038/s44294-026-00163-6](https://doi.org/10.1038/s44294-026-00163-6)

\* Equal contribution. <sup>†</sup> Joint senior authors. Open access, CC BY 4.0.

Cross-sectional study of 111 Swedish women (29 BV, 82 controls) combining shotgun
metagenomics, targeted LC-MS metabolomics and *Gardnerella* metagenome-assembled
genome (MAG) analysis, integrated with DIABLO and MOFA2.

## Data availability

| Data                                | Location                                                                |
| ----------------------------------- | ----------------------------------------------------------------------- |
| Raw sequencing reads                | ENA[PRJEB108308](https://www.ebi.ac.uk/ena/browser/view/PRJEB108308)     |
| Supplementary Data 1-4              | Zenodo[10.5281/zenodo.20734655](https://doi.org/10.5281/zenodo.20734655) |
| Supplementary Tables S1-S13         | `supplementary_tables/`                                               |
| Figures 1-7, S1                     | `figures/`                                                            |
| Processed omics matrices            | `objects/multiomics_views.rda`                                        |
| Sample group labels                 | `metadata/sample_metadata_public.rda`                                 |
| Individual-level clinical variables | Controlled access (see below)                                           |

**Zenodo contents.** Supplementary Data 1, the targeted LC-MS panel, is the assay
*method*: 191 compounds with retention times, MRM transitions and database
identifiers, one row per compound and no sample dimension. Supplementary Data 2
lists screened-but-excluded compounds, Data 3 the *Gardnerella* MAG quality
metrics, and Data 4 the full gene-metabolite correlation matrix (39 MB).

**Metabolomics levels.** The assay method is on Zenodo (Data 1). Processed
per-sample abundances, normalised, imputed and batch-corrected, are in
`objects/multiomics_views.rda` (111 x 133 positive mode, 111 x 31 negative mode).
Unprocessed per-sample intensities (`polar_positive_results.tsv`,
`polar_negative_results.tsv`, consumed by `00_data_processing_and_DIABLO.R`) are
held on institutional storage and are not deposited; contact the corresponding
author.

### Controlled access

Age, vaginal pH, antibiotic use and BV history are not distributed here. The
article reports them only as model covariates and aggregate results, and they are
keyed to the ENA sample identifiers, so releasing them would expose
per-participant health data beyond the scope of participant consent.

Requests go to the corresponding author, Kenny A. Rodriguez-Wallberg
(kenny.rodriguez-wallberg@ki.se), citing ethical approval 2021-05794-01. Access
requires a data transfer agreement and approval by the Swedish Ethical Review
Authority.

`metadata/sample_metadata_public.rda` pairs each sample identifier with its
BV/Control label. Release of these group labels is covered by ethical approval
2021-05794-01. The label carries no clinical measurement, is the study's primary
grouping variable, and is recomputable from the taxonomy matrix in
`objects/multiomics_views.rda` and from the public ENA metagenomes.

The de-identified objects shipped here are produced by
[`scripts/make_public_objects.R`](scripts/make_public_objects.R), which documents
what is removed and why.

## Repository structure

```
├── data/                                     Gardnerella pangenome (Panaroo, 61 MAGs)
├── metadata/sample_metadata_public.rda       meta_data: Sample + Group (n = 111)
├── objects/
│   ├── multiomics_views.rda                  X: metabolomics +/-, pathways, taxonomy
│   ├── group_labels.rda                      Y_aligned: BV / Control
│   ├── diablo_selected_features.rda          DIABLO feature selection
│   └── genome_metabolome_results_public.rda  genome-metabolome results
├── scripts/
│   ├── _metadata_access.R                    metadata loader and clinical gating
│   ├── 00_data_processing_and_DIABLO.R       upstream pipeline and DIABLO
│   ├── 01_MOFA2_integration.R                MOFA2 factor analysis
│   ├── 02_genome_metabolome.R                genome-metabolome integration, Figure 7
│   ├── make_public_objects.R                 de-identification (authors only)
│   └── qa/                                   privacy scanners
│       ├── scan_pii.sh                       static text scan
│       ├── scan_rda.R                        deep object scan
│       └── scan_values.R                     value-level scan (authors)
├── figures/                                  Figure1-7, FigureS1 (assembled; see below)
└── supplementary_tables/                     Table_S1-S13
```

## Supplementary tables

Numbering matches the published Supplementary Information.

| Table    | Contents                                                   | Regenerated by |
| -------- | ---------------------------------------------------------- | -------------- |
| S1 (a-d) | Targeted LC-MS panel and chromatography                    | curated        |
| S2       | Metabolite detection and missingness by BV status          | curated        |
| S3       | Antibiotic-exclusion sensitivity analysis                  | curated        |
| S4       | *Gardnerella* MAG quality and assembly metrics           | curated        |
| S5       | PERMANOVA of community composition                         | `00_*.R`     |
| S6       | ANCOM-BC2 differential abundance, taxa                     | `00_*.R`     |
| S7       | ANCOM-BC2 differential abundance, MetaCyc pathways         | `00_*.R`     |
| S8, S9   | limma differential metabolites, positive and negative mode | `00_*.R`     |
| S10      | DIABLO variance explained                                  | `00_*.R`     |
| S11      | DIABLO feature loadings                                    | `00_*.R`     |
| S12      | MOFA2 variance explained                                   | `01_*.R`     |
| S13      | MOFA2 complete factor loadings                             | `01_*.R`     |

`02_*.R` additionally writes `gene_metabolite_correlations_q0.1.xlsx`, which is
not a numbered supplementary table; the full matrix is Zenodo Supplementary
Data 4.

## Running the analyses

```
00_data_processing_and_DIABLO.R      (requires raw data, not included)
              │
      ┌───────┴───────┐
      ▼               ▼
01_MOFA2_integration  02_genome_metabolome
```

Run from the repository root so relative paths resolve:

```r
setwd("/path/to/mitch_manuscript")
source("scripts/01_MOFA2_integration.R")
source("scripts/02_genome_metabolome.R")
```

| Script                              | Purpose                                                                                        | Runs from this repository alone                |
| ----------------------------------- | ---------------------------------------------------------------------------------------------- | ---------------------------------------------- |
| `00_data_processing_and_DIABLO.R` | Read processing, MetaPhlAn4 profiling, decontam, ANCOM-BC2, LC-MS normalisation, limma, DIABLO | No; requires raw data from ENA and the authors |
| `01_MOFA2_integration.R`          | MOFA2 factor analysis across four views                                                        | Yes, except the age and pH sub-analyses        |
| `02_genome_metabolome.R`          | Mantel, Procrustes, db-RDA, gene-metabolite correlations, Figure 7                             | Yes                                            |

`scripts/_metadata_access.R` loads the participant-level metadata when present and
otherwise falls back to the public file. Analyses requiring the withheld variables
(age and pH summaries, the enhanced 4-view classification, pH pattern analysis,
factor-versus-pH correlation) are skipped with an explanatory message rather than
failing.

The public cache omits `dbrda_cond`, the db-RDA fitted with
`Condition(pH_z + Age_z + Antibiotic + BV_History)`, because its `$pCCA$QR`
component stores the QR decomposition of that conditioning matrix and would allow
the withheld values to be reconstructed. It ships the quantity the analysis uses,
`dbrda_cond_adjR2` = 0.0854, together with the permutation test `anova_cond`.

Of the published figures, only `Figure1.png` is written directly by a script
(`00_data_processing_and_DIABLO.R`). `02_genome_metabolome.R` writes the Figure 7
panels (`Figure7a_procrustes`, `Figure7b_varpart`, `Figure7c_heatmap` and a
composite), and `01_MOFA2_integration.R` writes MOFA2 panels; these generated
filenames are gitignored. The final `Figure2`-`Figure7` and `FigureS1` committed
here were assembled from those panels for publication.

Values reproduced from this repository: Mantel *r* = 0.148, *p* = 0.034;
Procrustes *r* = 0.385, *p* = 0.018; db-RDA conditioned adjusted *R*² = 0.085.

## Quality control

```bash
bash scripts/qa/scan_pii.sh      # static text scan
Rscript scripts/qa/scan_rda.R    # deep scan inside .rda, .rds and .xlsx
Rscript scripts/qa/scan_values.R # value-level scan (authors only)
```

`.gitignore` guards paths, not contents, so `scan_rda.R` recurses through lists,
S4 slots, attributes and environments, to a depth limit of 12, matching withheld
variable names against names, dimnames and factor levels.

Name matching cannot catch a renamed column, so `scan_values.R` compares the
actual withheld vectors against the numeric vectors in the publishable objects,
including individual matrix columns, matching exact values, sorted values and
perfect linear rescales. It needs the private inputs and therefore skips
automatically wherever they are absent, including in CI.

Neither scanner is a proof of absence. They are a regression guard against the
specific failure modes found while preparing this release: a clinical table
embedded in a results object, a conditioning matrix retained inside a model fit,
and values surviving under a renamed column or an attribute.

The first two run in `.github/workflows/privacy-check.yml`; all three run as a
pre-commit hook:

```bash
git config core.hooksPath .githooks
```

## Dependencies

**Bioconductor:** MOFA2, mixOmics, ANCOMBC, DESeq2, phyloseq, limma, decontam,
impute, sva, BiocParallel, ALDEx2, Maaslin2

**CRAN:** tidyverse, readxl, reshape2, data.table, vegan, ade4, compositions,
ggplot2, ggrepel, ggtext, pheatmap, corrplot, RColorBrewer, cowplot, patchwork,
gridExtra, scales, randomForest, caret, openxlsx, broom, matrixStats,
MetaboAnalystR

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("MOFA2", "mixOmics", "ANCOMBC", "DESeq2", "phyloseq",
                       "limma", "decontam", "impute", "sva", "BiocParallel",
                       "ALDEx2", "Maaslin2"))
install.packages(c("tidyverse", "readxl", "reshape2", "data.table", "vegan",
                   "ade4", "compositions", "ggrepel", "ggtext", "pheatmap",
                   "corrplot", "RColorBrewer", "cowplot", "patchwork",
                   "gridExtra", "randomForest", "caret", "openxlsx", "broom",
                   "matrixStats", "MetaboAnalystR"))
```

## Ethics

Approved by the Swedish Ethical Review Authority (Etikprövningsmyndigheten), Lund
Division 2 Medicine, approval 2021-05794-01, and conducted in accordance with the
Declaration of Helsinki. Written informed consent was obtained from all
participants before enrolment and sample collection.

## Citation

```bibtex
@article{Shabana2026MiTCH,
  author  = {Shabana, Hana and Dube, Faruk and Lahtinen, Emilia and
             Bergman R{\o}the, Emelie and Engstrand, Lars and
             Warchavchik Hugerth, Luisa and Schuppe-Koistinen, Ina and
             Rodriguez-Wallberg, Kenny A.},
  title   = {Multi-omics investigation of metabolomic and microbial features
             in bacterial vaginosis},
  journal = {npj Women's Health},
  volume  = {4},
  number  = {1},
  pages   = {26},
  year    = {2026},
  doi     = {10.1038/s44294-026-00163-6}
}
```

## License

Code in `scripts/` is MIT licensed. Figures, supplementary tables and processed
data objects are CC BY 4.0, matching the article. See [LICENSE](LICENSE).
Individual-level participant data are not covered by either licence and are not
distributed here.
