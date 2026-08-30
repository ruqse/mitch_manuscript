# Supplementary Table Legends

Column abbreviations used across tables: **log2FC**, log2 fold-change (BV vs control; positive = enriched
in BV); **fold_change**, linear fold-change (= 2^|log2FC|); **LFC**, log-fold change; **SE**, standard
error; **q / FDR / adj.P.Val / padj**, Benjamini–Hochberg false-discovery-rate-adjusted p-value;
**pval / P.Value / Fisher_p**, nominal p-value; **W_stat**, ANCOM-BC2 test statistic; **AveExpr**, average
(log2-scale) abundance across samples; **t**, moderated t-statistic (limma); **B**, log-odds of
differential abundance (limma); **ANI**, average nucleotide identity; **CDS**, protein-coding sequences;
**BER**, balanced error rate.

---


**Table S1. Targeted LC-MS metabolomics panel and chromatography.**
(a) Chromatographic gradient for positive-ion mode; (b) chromatographic gradient for negative-ion mode;
(c) retention time and MS/MRM parameters for each of the 191 reported compounds, with identification
level, MRM transitions, and database identifiers (KEGG, HMDB, PubChem, InChI, CAS, molecular formula and
weight); (d) the 72 polar compounds (46 positive, 26 negative ionization mode) screened but not included
among the final reported analytes, with the reason for exclusion (e.g., not detected, not reliably
quantifiable), reported with the same MS/transition and database-identifier fields as (c).

**Table S2. Metabolite detection and missingness by BV status.** Per-feature missingness on the raw
(pre-imputation) panel for the 111-cohort. Columns: Mode (pos/neg); Metabolite; Super_Pathway;
Chemical_Class; Pct_missing_overall; Pct_missing_BV; Pct_missing_Control; Fisher_p (Fisher exact test of
group-dependent detection); Fisher_FDR (BH-adjusted).

**Table S3. Antibiotic-exclusion sensitivity analysis.** Per-feature effect sizes from the main model
versus a model excluding the 11 participants with current antibiotic use. Sheet *Table S3a* (taxa): taxon;
log2FC_main; q_main; log2FC_noABX; q_noABX; same_direction. Sheet *Table S3b* (metabolites): mode (pos/neg);
Metabolite; log2FC_main; FDR_main; log2FC_noABX; FDR_noABX; same_direction. "_noABX" = antibiotic-excluded
cohort (n = 100). The taxonomic main column reproduces Table S6 exactly; the antibiotic-excluded column was
obtained by re-fitting ANCOM-BC2 (v2.4.0) on the same filtered taxon table and settings, dropping the
now-constant antibiotic term, so both columns share the same baseline and are directly comparable. log2FC_noABX is NA for taxa not testable in the reduced cohort (e.g. *L. portuensis*).

**Table S4. Gardnerella MAG quality and assembly metrics.** All 61 recovered Gardnerella MAGs. Columns:
MAG_ID; Sample (source sample); In_cohort (Yes = from a cohort participant, n = 35; No = n = 26);
BV_status (BV / Control / non-cohort); GTDB_Tk_species; closest_genome_ANI (%); binner; Completeness_pct
and Contamination_pct (CheckM2); Genome_Size_Mb; Contig_N50; Total_Contigs; Gene_count_CDS;
GC_pct. MAGs were retained at ≥70% completeness and ≤5% contamination.

**Table S5. PERMANOVA of community composition.** Partitioning of Bray–Curtis dissimilarity (adonis2 function, vegan R package).
Columns: Term (model term); Df (degrees of freedom); SumOfSqs (sum of squares); R2 (% variance explained);
F_stat (pseudo-F); P_value (999 permutations).

**Table S6. ANCOM-BC2 differential abundance of taxa (BV vs control).** Adjusted model
(~ Group + vaginal pH + antibiotic + BV history + age). Columns: taxon; log2FC_GroupBV (log2 fold-change
for the BV contrast, converted from ANCOM-BC2's natural-log LFC); fold_change (linear, 2^|log2FC|);
q_GroupBV (BH-adjusted p-value); passed_sensitivity (TRUE if retained in ANCOM-BC2's internal
sensitivity analysis; of the seven taxa at q < 0.05, five are also sensitivity-passing).

**Table S7. ANCOM-BC2 differential abundance of MetaCyc pathways (BV vs control).** Columns: pathway;
log2FC; fold_change; SE (log2 scale); W_stat; pval; padj; is_DA (TRUE if padj < 0.05 and |log2FC| ≥ 1);
passed_sensitivity (TRUE if retained in ANCOM-BC2's internal sensitivity analysis).

**Table S8 / Table S9. limma differential abundance of metabolites, positive (S8) and negative (S9) ion
mode (BV vs control).** Covariate- and batch-adjusted model. Columns: Metabolite; log2FC; fold_change
(2^|log2FC|); AveExpr; t; P.Value; adj.P.Val; B; FDR; direction (Higher in BV / Higher in Control).
Fold-changes are derived from the un-scaled natural-log limma coefficient; significance is unchanged by
this rescaling.

**Table S10. DIABLO variance explained.** Columns: Block (omics block); Component (1–3);
Variance_Explained (% of block variance captured by the component).

**Table S11. DIABLO feature loadings.** Columns: Feature; Loading (signed component loading); Component;
Block; Abs_Loading (absolute loading, used for ranking).

**Table S12. MOFA2 variance explained.** Percentage of variance explained by each MOFA2 factor (columns)
within each omics view (rows).

**Table S13. MOFA2 complete factor loadings.** Per-feature loadings for all factors. Columns: View;
Feature; Factor1–Factor7; Max_Abs_Loading; DIABLO_Overlap (whether the feature was also selected by
DIABLO).

**Table S14. Gene–metabolite correlations (Gardnerella shell genome).** Spearman correlations between
shell-genome gene clusters (15–95% prevalence) and metabolites. Columns: Gene_ID; Metabolite; rho
(Spearman coefficient); pvalue; Gene_Name; Gene_Annotation; qvalue (BH-adjusted).
