# GENESIS_iClust
Single-cell iterative clustering pipeline for GENESIS

Forked from [GENESIS_iClust](https://github.com/hoondy/GENESIS_iClust) v1.0.1 (commit `638116a`) and updated to v1.3.0: v1.1.0 and v1.2.0 (commit `727f1e8`) change step 10, and v1.3.0 (commit `42ad62d`) changes gene names in step 02. The pipeline version is `__version__` in `00_shared.py`.

The statistical methods, design and software versions are described in [METHODS.md](METHODS.md). This README covers how to run the code.

## Steps

| Script | Step |
| --- | --- |
| `01_aggregate.py` | Merge per-sample h5ad files on disk |
| `02_anno.py` | Add Ensembl gene annotation, index genes by symbol (one name per gene) |
| `03_qc.py` | Gene, cell and channel QC |
| `04_norm.py` | Shifted log1p normalisation (raw counts kept in `raw.X`) |
| `05_score.py` | Gene-signature scores (cell cycle, gender, ...) |
| `06_iclust.py` | Iterative Leiden clustering: class > subclass > subtype |
| `07_scrublet.py` | Doublet detection and removal |
| `08_subset.py` | Remove listed clusters |
| `09_annotate.py` | Apply curated annotation and re-embed |
| `10_clean.py` | Remove lineage-mixed cells, and cells that sit among another subclass |
| `11_umap.py` | Final embedding and UMAPs |

The Snakemake files chain the steps:
- `TEST_pass2.smk`: 01-07, then a second 06 pass
- `TEST_pass3.smk`: 08 and a third 06 pass
- `TEST_rc2.smk`: 09-11

## Gene names

Step 02 names each gene using `annotation/all_genes_EnsemblInfo_104_1to1.csv`. In this table every Ensembl ID has its own name: its HGNC symbol, or its Ensembl ID if it has no symbol or shares its symbol with another gene (`02_anno.py` says which gene keeps a shared symbol). If two genes would still end up with the same name, step 02 stops with an error. Copy the table into each cohort's `annotation/` folder, where `TEST_pass2.smk` looks for it.

Up to rc2, step 02 used the full Ensembl 104 table (`annotation/all_genes_EnsemblInfo_104.csv`). That table gives a few symbols to two genes, and step 02 renamed the second gene with `#1`, so a symbol could point to the wrong gene: in rc2, `TBCE` is a readthrough gene and the real TBCE is `TBCE#1`. Eight genes are therefore named differently from rc3 onward. To compare releases, match genes on `gene_id`; the table's `rc2_var_name` column gives each gene's rc2 name.

| Ensembl ID | rc2 | rc3 onward |
| --- | --- | --- |
| ENSG00000284770 | `TBCE#1` | `TBCE` |
| ENSG00000285053 | `TBCE` | `ENSG00000285053` |
| ENSG00000145075 | `CCDC39#1` | `ENSG00000145075` |
| ENSG00000206195 | `DUXAP8#1` | `DUXAP8` |
| ENSG00000261186 | `LINC01238#1` | `LINC01238` |
| ENSG00000280709 | `LINC02203#1` | `LINC02203` |
| ENSG00000258724 | `PINX1#1` | `PINX1` |
| ENSG00000268500 | `SIGLEC5#1` | `SIGLEC5` |

## Shared code

- [SCTools](https://github.com/hoondy/SCTools) provides the general-purpose helpers: on-disk h5ad I/O, HVF selection, QC thresholds, signature scores, outlier and lineage-mixing detection. Install SCTools ≥ 0.8.0 in the pipeline environment (`py312`); HVF selection uses the `min_cells` option of `sct.pp.scanpy_hvf_h5ad`, and step 10 the `score="nuclei"` option of `sct.tl.detect_lineage_mixing` and the `method="neighbor_label"` option of `sct.tl.detect_outliers`.
- HVFs (06, 09, 11) are always picked with Scanpy's `cell_ranger` flavor from protein-coding autosomal genes that are expressed in the input. Genes with no counts are dropped first, because `cell_ranger` fails with "Bin edges must be unique" on small subclusters when many genes are unexpressed.
- `00_shared.py` holds the pipeline-specific helpers that several steps reuse: shared defaults, HVF loading, the PCA/Harmony/kNN/UMAP embedding, and UMAP plotting with the GENESIS palette. Steps load it with `importlib.import_module('00_shared')`, so it must sit in the same directory as the numbered scripts.

## Releases and pipeline versions

Step 10 and gene names changed after rc2, so rc2 cannot be reproduced with the current defaults.

| Release | Pipeline | Step 02: gene names | Step 10: mixed cells | Step 10: outliers |
| --- | --- | --- | --- | --- |
| rc2 | v1.0.1 (SCTools 0.5.0) | full Ensembl 104 table; a shared symbol's second gene renamed `#1` | share of marker RNA from other classes ≥ 0.3 | local outlier factor (LOF) per subclass |
| rc3 onward | v1.3.0 | one name per gene (v1.3.0) | RNA from any one other class ≥ 0.25 of a nucleus (v1.1.0) | at most 4 of the 15 nearest cells share the subclass (v1.2.0) |

LOF compared each cell only with its own subclass. It missed cells that sit inside another subclass's island on the UMAP, which were the outliers it was meant to remove, and it trimmed small subtypes and donor-specific cell states instead. `--outlier_method lof` still runs the rc2 outlier rule.

