# MitCH Manuscript: Multi-omics Analysis of Bacterial Vaginosis

## Study

**Bacterial vaginosis as a metabolic state: multi-omics integration reveals community-driven polyamine signatures independent of *Gardnerella* strain diversity**

Hana Shabana\*, Faruk Dube\*, Emilia Lahtinen, Emelie Bergman Rothe, Lars Engstrand, Luisa W. Hugerth, Ina Schuppe-Koistinen#, Kenny A. Rodriguez-Wallberg#

\*Equal contribution as first authors; #Equal contribution as senior authors

## Description

Companion code and processed data for the MitCH multi-omics BV study. This repository contains the analysis scripts, processed intermediate objects, and supplementary materials for a cross-sectional study of 111 Swedish women (29 BV-positive, 82 controls) integrating shotgun metagenomics, targeted LC-MS metabolomics, and *Gardnerella* pangenomics.

## Data Availability

| Data type | Location | Status |
|-----------|----------|--------|
| Raw sequencing reads | ENA accession [PRJEB108308](https://www.ebi.ac.uk/ena/browser/view/PRJEB108308) | Available upon publication |
| Processed objects and metadata | This repository | Available |
| Large supplementary files | Zenodo (DOI to be added) | Available upon publication |
| Raw metabolomics data | Available from authors upon request | - |

## Repository Structure

```
mitch_manuscript/
├── README.md
├── data/                           # Gardnerella pangenome gene presence/absence
├── metadata/
│   └── sample_metadata.rda         # Sample metadata (n=111)
├── objects/                        # Cached R objects for reproducibility
│   ├── multiomics_views.rda        # 4-view omics data (metabolomics+/-, pathways, taxonomy)
│   ├── group_labels.rda            # BV/Control group labels
│   ├── diablo_selected_features.rda # DIABLO feature selection results
│   └── genome_metabolome_results.rda # Genome-metabolome integration results
├── scripts/
│   ├── 00_data_processing_and_DIABLO.R  # Raw data processing and DIABLO integration
│   ├── 01_MOFA2_integration.R           # MOFA2 multi-omics factor analysis
│   └── 02_genome_metabolome.R           # Genome-metabolome integration + Figure 7
├── figures/                        # Manuscript figures
└── supplementary_tables/           # Tables S1-S11
```

## Script Execution

### Execution order

```
00_data_processing_and_DIABLO.R  (requires raw data, not included)
         |
         v
   +-----+-----+
   |             |
   v             v
01_MOFA2       02_genome_metabolome
(independent)  (independent)
```

### Script descriptions

| Script | Purpose | Raw data needed? |
|--------|---------|-----------------|
| `00_data_processing_and_DIABLO.R` | Processes raw sequencing and metabolomics data, performs differential abundance (ANCOM-BC2, limma), and fits DIABLO multi-omics integration | Yes (see script header for data sources) |
| `01_MOFA2_integration.R` | Unsupervised multi-omics factor analysis (MOFA2) with 4 data views | No (uses cached objects) |
| `02_genome_metabolome.R` | *Gardnerella* pangenome-metabolome integration (Mantel, Procrustes, db-RDA) and publication Figure 7 | No (uses cached objects) |

**Reproducibility note:** Scripts 01 and 02 are fully reproducible using the cached objects in `objects/` and `metadata/`. Script 00 documents the complete upstream pipeline but requires raw data files available from ENA and from the authors.

## R Dependencies

### Bioconductor packages
MOFA2, mixOmics, ANCOMBC, DESeq2, phyloseq, limma, decontam, impute, sva, BiocParallel, ALDEx2, Maaslin2

### CRAN packages
tidyverse, dplyr, tidyr, purrr, readr, tibble, magrittr, stringr, readxl, reshape2, data.table, ggplot2, ggrepel, ggtext, pheatmap, corrplot, RColorBrewer, cowplot, patchwork, gridExtra, scales, vegan, ade4, compositions, randomForest, caret, openxlsx, broom, matrixStats, MetaboAnalystR

### Installation

```r
# Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("MOFA2", "mixOmics", "ANCOMBC", "DESeq2", "phyloseq",
                        "limma", "decontam", "impute", "sva", "BiocParallel",
                        "ALDEx2", "Maaslin2"))

# CRAN
install.packages(c("tidyverse", "readxl", "reshape2", "data.table", "vegan",
                    "ade4", "compositions", "ggrepel", "ggtext", "pheatmap",
                    "corrplot", "RColorBrewer", "cowplot", "patchwork",
                    "gridExtra", "randomForest", "caret", "openxlsx",
                    "broom", "matrixStats", "MetaboAnalystR"))
```

## Ethics

This study was approved by the Swedish Ethical Review Authority (reference: 2021-05794-01). All participants provided written informed consent.

## Citation

*Manuscript under review*

## License

MIT
