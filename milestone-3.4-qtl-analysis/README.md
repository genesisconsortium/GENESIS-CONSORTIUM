# Milestone 3.4: QTL analysis

## Freeze 1 pseudobulk preprocessing

This page documents the preprocessing steps used to generate the GENESIS Freeze 1 pseudobulk Pearson residual matrices for downstream sc-eQTL, sc-eTR, pan-region QTL, region-specific QTL, and TWAS analyses.

Freeze 1 currently includes:

- 3,139 final WGS-ready donor-region samples
- 2,420 unique primary WGS genotypes
- 16,872,844 nuclei
- 7 brain regions

## 1. Input h5ad files

h5ad inputs are read from:

`/sc/arion/projects/CommonMind/genesis/pseudobulk_dreamlet_rc2/config/pseudobulk_inputs.csv`

Each row in this config file defines one dataset-region analysis, including project, cohort, h5ad path, donor column, region column, default region, metadata cohort, and output directory.


## 2. Sample metadata and primary genotype assignment

Metadata file:

`/sc/arion/projects/CommonMind/genesis/metadata/outputs/samples_single_cell_primary_GT.csv`

Final two-column sample-to-WGS mapping file:

`/sc/arion/projects/CommonMind/genesis/metadata/primary_gt_for_scRNAseq.csv`

Notes:

- `individualID` was used as the h5ad donor ID.
- `primary_genotype` was used as the final WGS genotype ID.
- Samples without a primary genotype were removed.
- When multiple WGS files were available for the same donor, the primary genotype was selected using the highest available WGS coverage/read-depth metric from Picard technical covariates.

## 3. Pseudobulk aggregation

Counts were aggregated to donor-level pseudobulk separately by brain region and annotation level.

Annotation levels generated:

- `class`
- `subclass`

## 4. Covariates and residualization

Pseudobulk expression was residualized using the model:

```r
~ sex + scale(age) + log(n_genes) + percent_mito +
  mito_genes + mito_ribo + ribo_genes
```

## 5. Outlier QC

Class-level pseudobulk QC was used for donor-region outlier detection.

Outlier detection steps:
- residuals from step 4 were rank-normalized
- PCA was performed on residualized matrices
- donor-region outlier testing was run in PCA space
- samples with FDR < 1e-5 were removed
- class-level outliers were removed from both class and subclass outputs

Five donor-region outliers were removed before final residual generation.

