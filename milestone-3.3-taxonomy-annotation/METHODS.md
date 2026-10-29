# Methods: GENESIS unified single-cell RNA-seq processing pipeline

This document describes the statistical methods, design and processing steps of the GENESIS unified processing pipeline, version 1.0.1 (the code in this directory, forked from GENESIS_iClust v1.0.1). Instructions for running the code are in [README.md](README.md). Parameter values given below are the defaults used by the Snakemake workflows in this directory. Throughout, "cell" refers to one barcode, whether it captured a whole cell or a nucleus.

## Pipeline design

The pipeline turns the per-library count matrices of one cohort into a quality-controlled, normalised, batch-integrated and hierarchically annotated dataset. It consists of eleven Python steps (`01`–`11`) run by three Snakemake workflows (Table 1). Its design follows six principles:

1. **Cohort-level processing.** All libraries of a cohort are merged before quality control (QC), so QC thresholds, normalisation and feature statistics are estimated on the whole cohort.
2. **Iterative clustering.** Cells are clustered at three levels (class > subclass > subtype). Each population is re-analysed from feature selection onwards, so variation within a population drives its embedding instead of differences between populations.
3. **Three clustering passes.** Passes are separated by doublet removal and by curated exclusion of clusters that do not replicate in independent datasets, so the final clusters are derived from cleaned data.
4. **One procedure at every level.** Feature selection, dimensionality reduction, covariate regression, batch integration, graph construction and visualisation are the same at every level and pass (implemented once in `00_shared.py`).
5. **On-disk operations.** Matrices are stored as compressed sparse row (CSR) data in H5AD files. Concatenation and subsetting stream rows between files in chunks of 500,000 cells, so memory use does not scale with cohort size.
6. **Fixed random seeds.** Every stochastic step uses random seed 0.

**Table 1 | Processing steps.**

| Workflow | Step | Script | Output |
| --- | --- | --- | --- |
| `TEST_pass2.smk` | Aggregation | `01_aggregate.py` | merged count matrix |
| | Gene annotation | `02_anno.py` | genes indexed by symbol, with Ensembl annotation |
| | Quality control | `03_qc.py` | QC-filtered cells and genes |
| | Normalisation | `04_norm.py` | log-normalised matrix, raw counts kept |
| | Signature scores | `05_score.py` | cell-cycle and other signature scores |
| | Clustering, pass 1 | `06_iclust.py` | class/subclass/subtype labels |
| | Doublet removal | `07_scrublet.py` | singlets |
| | Clustering, pass 2 | `06_iclust.py` | labels on singlets |
| `TEST_pass3.smk` | Cluster exclusion | `08_subset.py` | listed subtypes removed |
| | Clustering, pass 3 | `06_iclust.py` | final subtypes |
| `TEST_rc2.smk` | Taxonomy assignment | `09_annotate.py` | curated class/subclass labels, re-embedded |
| | Mixed-lineage and outlier removal | `10_clean.py` | cleaned cells |
| | Final embedding | `11_umap.py` | release dataset |

## Input data

The input is one H5AD file per sequencing library. Each file holds a cell × gene matrix of raw unique molecular identifier (UMI) counts, with genes indexed by Ensembl gene ID. All libraries of a cohort must share the same gene set in the same order, i.e. they must come from one reference annotation. The cell metadata must contain the library identifier (`libraryID`) or, failing that, the donor identifier (`individualID`). It must also contain the variable used as the batch for integration, which is `libraryID` in the GENESIS workflows. Upstream processing (read alignment, UMI counting and sample demultiplexing) is outside the scope of this pipeline.

## Aggregation and gene annotation

The libraries listed in a manifest were concatenated row-wise on disk (`sct.io.concat_on_disk`). The metadata of all libraries were concatenated, and the count matrix of each library was appended to a single CSR matrix without being loaded into memory as a whole.

Genes were annotated from the Ensembl release 104 gene table for GRCh38 (`annotation/all_genes_EnsemblInfo_104.csv`). The table supplied, for each Ensembl gene ID, the HGNC symbol, gene biotype, chromosome, and start and end position. Genes without an HGNC symbol kept their Ensembl ID as their name. Genes were then indexed by symbol, and duplicate symbols were made unique by appending `#1`, `#2`, … to later occurrences.

## Quality control

**Gene filtering.** QC was performed jointly on all cells of the cohort (pegasus `identify_robust_genes`).
- Genes detected (≥1 UMI) in no cell were removed.
- Of the remaining genes, those detected in at least 0.05% of cells ("robust" genes) were retained for all downstream analyses.
- Retained genes were flagged as:
  - ribosomal (symbol beginning with `RP`)
  - mitochondrially encoded (symbol beginning with `MT-`)
  - protein coding (Ensembl biotype `protein_coding`)
  - nuclear-encoded mitochondrial, if listed in MitoCarta3.0 (ref. 1)
- Two further flags were derived:
  - `robust_protein_coding`: robust, protein-coding genes that are not ribosomal, mitochondrially encoded or MitoCarta3.0 genes
  - `robust_protein_coding_autosome`: the same genes restricted to chromosomes 1–22

These flags were stored for downstream use; they were not used to filter genes.

**Cell filtering.** Three metrics were computed per cell over the retained genes: the total UMI count (`n_counts`), the number of detected genes (`n_genes`), and the percentage of UMIs from mitochondrially encoded genes (`percent_mito`). (pegasus keeps `n_counts` and `n_genes` unchanged if they are already present in the input.)

Thresholds for `n_counts` and `n_genes` were set by a median absolute deviation (MAD) rule on the log scale (`sct.pp.qc_boundary`). For a metric $x$, let $y = \log(1 + x)$, $m = \operatorname{median}(y)$ and $\mathrm{MAD} = \operatorname{median}(|y - m|)$, unscaled. The lower and upper thresholds were then $L = \exp(m - k \cdot \mathrm{MAD})$ and $U = \exp(m + k \cdot \mathrm{MAD})$, with $k = 3$.

A cell passed QC if all of the following held:
- $L \le$ `n_counts` $< U$
- $L \le$ `n_genes` $< U$ (with that metric's own thresholds)
- `percent_mito` $< 5$

**Library filtering.** Libraries (`libraryID`, or `individualID` when the library identifier was absent) with fewer than 50 cells passing QC were removed entirely. Cells failing any criterion were removed, and the per-step filter statistics and diagnostic histograms were saved.

## Normalisation

Raw counts were kept as `raw.X`. Expression was normalised by shifted logarithm (scanpy `normalize_total` with `target_sum=None`, then `log1p`; ref. 2):

$$\tilde{x}_{ig} = \log\!\left(1 + \frac{x_{ig}\, M}{n_i}\right),$$

where $x_{ig}$ is the UMI count of gene $g$ in cell $i$, $n_i$ is the total count of cell $i$, and $M$ is the median of $n_i$ across all cells. Logarithms are natural.

## Gene-signature scores

The QC metrics `n_counts`, `n_genes` and `percent_mito` were not recomputed after normalisation. They remain the raw-count values from the QC step, and these are the values used as covariates below.

Gene-signature scores were computed with pegasus `calc_signature_score` for its predefined human signatures:
- cell cycle: `G1/S` and `G2/M` genes from Tirosh et al. (ref. 4), with symbols updated as in Seurat's `cc.genes.updated.2019`
- sex: `female_score`, `male_score`
- mitochondrial: `mito_genes`, the 13 protein-coding chrM genes; and `mito_ribo`, nuclear-encoded mitochondrial ribosomal genes
- ribosomal: `ribo_genes`
- apoptosis: `apoptosis`, from KEGG

The scoring is an extension of the method of Jerby-Arnon et al. (ref. 5):
1. Genes were assigned to 50 equal-frequency bins of mean expression.
2. Each gene's log-normalised expression was centred on its mean and on a cell-specific background expectation, estimated from all genes in the same expression bin.
3. The centred value was scaled by the cell-specific background dispersion of that bin.
4. A signature score is the average of these standardised values over the signature's genes.

Derived quantities:
- `cycle_diff` = `G2/M` − `G1/S` (used as a covariate below)
- `cycling` = max(`G1/S`, `G2/M`)
- `predicted_phase`: cells in the higher of two *k*-means clusters of `cycling` were assigned G1/S or G2/M by their larger score; the remaining cells were labelled G0
- `gender_score` = `male_score` − `female_score`
- `predicted_gender`: female / uncertain / male, from three *k*-means clusters of `gender_score`

## Feature selection

Highly variable features (HVFs) were selected anew for every clustering and embedding run (`00_shared.read_input_with_hvf`, which calls `sct.pp.scanpy_hvf_h5ad`).

**Candidate genes.** Candidates were protein-coding genes outside chromosomes X, Y and MT that were detected in at least one cell of the input being analysed. Genes without counts in the input were removed first (`min_cells=1`). In small subsets, such genes can make up more than 15% of genes; their identical zero means then produce duplicate bin edges in the cell_ranger procedure, and scanpy raises "Bin edges must be unique". Genes without counts cannot be variable, so removing them changes no selection that would otherwise succeed.

**Selection.** HVFs were selected by the dispersion-based `cell_ranger` procedure (ref. 6) in scanpy:
1. For each gene, the mean $\mu_g$ and variance $\sigma_g^2$ of log-normalised expression were computed across the input's cells, giving dispersion $d_g = \sigma_g^2 / \mu_g$.
2. Genes were binned by $\mu_g$ at the 10th, 15th, …, 100th percentiles. The lowest bin holds 10% of genes and each further bin 5%.
3. Each dispersion was normalised within its bin as $(d_g - \operatorname{median}_b d)/\operatorname{MAD}_b(d)$, where the MAD is scaled for consistency with the normal distribution.
4. The genes with the highest normalised dispersion were selected; genes tied at the cut-off were also kept.

**Number of HVFs.**
- iterative clustering: 3,000 at class level, 2,000 at subclass level, 1,000 at subtype level
- re-embedding of annotated data (steps 09 and 11): 6,000

## Dimensionality reduction, covariate regression and batch integration

All operations below used pegasus (ref. 3) and were applied separately to every input (`00_shared.run_embedding`):

1. **Standardisation and PCA.** The log-normalised expression of the HVFs was standardised to zero mean and unit variance per gene, using that input's statistics. Values were clipped to [−10, 10]. Principal component analysis (PCA; scikit-learn, ref. 7) retained 30 components.
2. **Covariate regression.** Each principal component was regressed by ordinary least squares, with an intercept, on three per-cell covariates: the raw-count `n_counts` and `percent_mito` from QC, and `cycle_diff` from signature scoring. The residuals were kept as the corrected components.
3. **Batch integration.** The corrected components were integrated across batches with Harmony (ref. 8; harmony-pytorch implementation). The batch variable was `libraryID`, i.e. each sequencing library was its own batch. Settings:
   - diversity penalty θ = 2, soft-clustering width σ = 0.1, ridge penalty λ = 1, block proportion 0.05
   - number of Harmony clusters: min(100, ⌊*n*/30⌋) for *n* cells
   - at most 20 Harmony iterations (tolerance 10⁻⁴); at most 200 clustering iterations per round (tolerance 10⁻⁵)

The 30 Harmony-corrected components (`X_pca_regressed_harmony`) were the input to all neighbourhood-based steps.

## Neighbour graph, clustering and visualisation

**Neighbour graph.** A *k*-nearest-neighbour (kNN) graph was built on the Harmony-corrected components with Euclidean distance.
- *K* = min(100, ⌊√*n*⌋), counting the cell itself, so each cell had *K* − 1 neighbours. The cap therefore applies to inputs of fewer than 10,000 cells.
- Neighbours were found approximately with hnswlib (ref. 9; *M* = 20, ef_construction = 200, ef = 200). For inputs of at most 1,000 cells, exact neighbours were computed with scikit-learn.

**Edge weights.** Edges were weighted with a locally scaled Gaussian kernel (ref. 10):

$$w_{ij} = \sqrt{\frac{2\sigma_i\sigma_j}{\sigma_i^2 + \sigma_j^2}}\, \exp\!\left(-\frac{d_{ij}^2}{\sigma_i^2 + \sigma_j^2}\right),$$

where $d_{ij}$ is the distance between cells $i$ and $j$ and $\sigma_i$ is the median distance from cell $i$ to its neighbours. The graph was then:
1. symmetrised (union of directed edges)
2. density-normalised as $w_{ij}/(z_i z_j)$, where $z_i = \sum_j w_{ij}$ (ref. 11)
3. scaled to a median edge weight of 1 and rounded to two decimals; edges of weight 0 were dropped

**Clustering.** Clusters were found with the Leiden algorithm (ref. 12):
- leidenalg, with the `RBConfigurationVertexPartition` quality function on the weighted graph
- iterated until the partition no longer improved
- resolution γ = 0.1, 0.2 and 0.3 at the class, subclass and subtype levels

Clusters were numbered from 1 in order of decreasing size.

**Visualisation.** Two-dimensional embeddings were computed with UMAP (ref. 13; umap-learn):
- 15 neighbours, minimum distance 0.5, spread 1, spectral initialisation, on the Harmony-corrected components
- for inputs of at least 4,096 cells, the 15-nearest-neighbour graph came from the hnswlib index above; for smaller inputs, umap-learn computed exact neighbours itself

## Iterative clustering

Each clustering pass (`06_iclust.py`) labelled cells at three levels:

1. **Class.** All cells were processed as above (3,000 HVFs, γ = 0.1). This gave class labels `1`, `2`, ….
2. **Subclass.** For each class, its cells were written to a new file on disk and re-processed from feature selection onwards (2,000 HVFs, γ = 0.2). Subclass labels were prefixed with the class label: `<class>_1`, `<class>_2`, ….
3. **Subtype.** Each subclass was re-processed the same way (1,000 HVFs, γ = 0.3). This gave subtype labels `<subclass>_1`, ….

At every level, HVFs, standardisation, principal components, regression, Harmony integration, neighbour graph, clustering and UMAP were all recomputed within the subset. The labels of all three levels were written back to the full dataset, and the embedding of the class level was kept. Because PCA retains 30 components, every population that is re-clustered must contain at least 30 cells.

The pipeline ran three passes:
- **Pass 1:** on all QC-filtered cells
- **Pass 2:** on singlets, after doublet removal
- **Pass 3:** after curated subtypes were excluded

## Doublet detection

Doublets were detected on the output of pass 1 with the Scrublet algorithm (ref. 14), as implemented in pegasus `infer_doublets`. Scrublet was run independently for each library; libraries with fewer than 100 cells were not scored. Within each library:

1. **Features.** Raw counts were taken for the 3,000 HVFs selected at the class level of pass 1.
2. **Simulated doublets.** Twice as many artificial doublets as observed cells were simulated by summing the count profiles of randomly chosen pairs of observed cells.
3. **Embedding.** Observed and simulated profiles were normalised to their total counts and standardised using the observed cells' statistics. Thirty principal components were fitted on the observed cells, and the simulated doublets were projected onto them.
4. **Neighbours.** For each profile, the $k_{\mathrm{adj}} = \operatorname{round}(k(1 + r))$ nearest neighbours were found, where $k = \operatorname{round}(0.5\sqrt{n})$ for *n* observed cells and $r = 2$ is the simulation ratio.
5. **Doublet score.** Let $k_d$ be the number of simulated doublets among a profile's neighbours, $q = (k_d + 1)/(k_{\mathrm{adj}} + 2)$, and $\rho$ the expected doublet rate. The score was $\frac{q\rho/r}{(1-\rho) - q(1 - \rho - \rho/r)}$.
6. **Expected doublet rate.** $\rho$ came from a Poisson model of the 10x Genomics multiplet table. With $\lambda = 0.00785 \times n/500$, $\rho = \frac{1 - (1+\lambda)e^{-\lambda}}{1 - e^{-\lambda}}$.
7. **Threshold.** The doublet-score threshold was inferred automatically for each library from the kernel density of log-transformed simulated-doublet scores, using a signed-curvature criterion. Cells scoring above it were predicted doublets.

**Cluster-level test.** Predicted doublets were tested for enrichment in each pass-1 subtype with Fisher's exact test, corrected for multiple testing by the Benjamini–Hochberg procedure (FDR 0.05; ref. 15). Subtypes with a significantly higher doublet fraction than the overall rate were reported to inform curation; they were not removed automatically. All predicted doublets were removed, and the remaining singlets were re-clustered in pass 2.

## Curation: cluster exclusion and taxonomy assignment

Clusters were curated by testing whether each subtype replicates in independent datasets with unsupervised MetaNeighbor (refs. 16,17). The curation code is in the sibling directory [`milestone-3.3-taxonomy-integration`](../milestone-3.3-taxonomy-integration). It was run interactively on CAVATICA Data Studio and is not part of the Snakemake workflows. Its results were used to write the two curation inputs of this pipeline: the list of pass-2 subtypes to exclude, and the table assigning each pass-3 subtype a class and subclass.

**Replicability scores.** In each comparison, the GENESIS subtypes (raw counts) and the cell types of the reference datasets were combined on their shared genes. MetaNeighbor (pyMN, fast implementation, with node-degree normalisation) then proceeded as follows:
1. **Variable genes.** Within each dataset, genes were divided into ten equal-frequency bins of median expression, and the genes in the top quartile of variance within their bin were selected. Only genes selected in every dataset were used.
2. **Cell similarity.** Each cell's expression of these genes was ranked, centred and scaled to unit norm, so that the dot product of two cells $i$ and $j$ is their Spearman correlation $\rho_{ij}$.
3. **Votes.** For a training cell type $A$ in dataset $D_A$ and a cell $i$ of another dataset, the vote was

   $$v_i(A) = \frac{\sum_{j \in A} (1 + \rho_{ij})}{\sum_{j \in D_A} (1 + \rho_{ij})}.$$

4. **All-vs-all AUROC.** This AUROC measured, within each test dataset, how well $v(A)$ separated the cells of test type $B$ from the other cells of that dataset. The score of a pair was the mean of the two directions ($A$ training and $B$ testing, and the reverse).
5. **One-vs-best AUROC.** For each training type and test dataset, the five test types with the highest all-vs-all AUROC were compared pairwise, using only their own cells. The type whose cells received higher votes than every contender was the best match. Its score $s$ was its AUROC against the closest contender, which received $1 - s$.
6. **Top hits.** The top hit of a cell type in another dataset was its best match there, provided its all-vs-all AUROC was at least 0.8.

**Stage 1: GENESIS reference cohorts.** Subtypes of three GENESIS cohorts were compared with a cortical reference of five datasets. The three cohorts and their passes were GEN_A1 (pass 2), GEN_A2 (pass 3) and GEN_A3 (pass 2). The reference datasets were:
- SEA-AD (ref. 18): a random quarter of the cells
- Ma et al. (ref. 19): all cells
- Siletti et al. (ref. 20): all cells except the superclusters "Eccentric medium spiny neuron", "Hippocampal CA1-3", "Fibroblast" and "Splatter"
- GEN_A1 (PsychAD): a random 1/22 of the cells of each of its four data files
- GEN_A3 (AMP PD): a random quarter of the cells from the cortical regions PFC, PMC and PVC

Each reference kept its own published or previously released subclass labels. The combined reference had 1,251,229 cells and 15,112 genes. Random sampling used seed 42.

For each subtype, its top hit in each reference was recorded, and subtypes with the same combination of top hits were collapsed into one pattern. A curator assigned each pattern either a subclass of the GENESIS taxonomy or a mark for removal. Typical reasons for removal were a pattern specific to one cohort, few cells, or matches to several subclasses. Microglia and perivascular macrophages were combined into one subclass (`Micro_PVM`). The result mapped every retained subtype of GEN_A1–A3 to a curated subclass.

**Stage 2: further cohorts.** Subtypes of each further GENESIS cohort were compared with the curated subclasses of three reference cohorts:
- GEN_A1: a random half of the cells
- GEN_A2: all cells
- GEN_A3: the cells from the regions PFC, PMC and PVC

Subtypes removed in Stage 1 were left out of the references. A subtype $t$ was assigned subclass $c$ by three rules:
1. **All-vs-all rule.** Take the top hit of $t$ in each reference cohort and order the hits by AUROC. $c$ is the subclass of the highest hit. Let $m$ be the number of reference cohorts that contain $c$. The first $\min(m, \text{number of hits})$ hits must all be $c$.
2. **One-vs-best rule.** Apply the same rule to the reference subclasses that $t$ matched with a one-vs-best AUROC above 0.6.
3. **Agreement.** Both rules must select the same subclass $c$.

Subtypes that satisfied all three rules were assigned $c$. Subtypes that failed, including those with no top hit, were flagged as not replicated. For cohorts curated after GEN_A3 was split by brain region, the references were the released subclasses of GEN_A1 and GEN_A2.

**Pipeline inputs.**
1. **Exclusion.** Subtypes that were marked for removal or not replicated were listed in the exclusion file. The cells of the listed subtypes were removed (`08_subset.py`), and the remaining cells were re-clustered in pass 3.
2. **Taxonomy assignment.** Each retained pass-3 subtype was entered in an annotation table (columns `class`, `subclass`, `subtype`). Its subclass was the one assigned above. Its class was the parent of that subclass in the GENESIS taxonomy (`parent` column of `annotation/260122_GENESIS_color_palette.csv`, which also sets the class and subclass colours). Cells of subtypes absent from the table were removed.
3. **Re-embedding.** Annotated cells were re-embedded without clustering (`09_annotate.py`), using the same procedure with 6,000 HVFs.

## Removal of mixed-lineage cells

Cells whose transcriptome is a mixture of cell types from different lineages were identified with a multinomial mixture model (`sct.tl.detect_lineage_mixing`; SCTools). Such mixtures arise from heterotypic doublets, nuclei carrying material from another cell, and heavy ambient contamination.

**Groups and lineages.**
- Mixture components were the curated subclasses; lineages were the curated classes.
- OPC and Oligo were treated as one lineage, so that genuine transitions between oligodendrocyte precursors and oligodendrocytes are not counted as mixing.
- Subclasses with fewer than 20 cells were not modelled, and their cells were not scored.

**Reference cells and marker panel.** All calculations used raw counts.
1. Up to 5,000 cells per subclass were sampled at random as reference cells.
2. For each subclass and each lineage, the mean of per-cell gene proportions over reference cells was computed.
3. A marker panel was assembled from up to 150 genes per subclass and up to 150 per lineage. A gene was eligible if:
   - its log2 ratio of mean proportion to the highest mean among other groups was at least 1, and
   - its mean proportion was at least the median, across genes, of the per-gene maximum subclass mean.
4. Mitochondrially encoded (`MT-`) and ribosomal-protein (`RPL`/`RPS`) genes were excluded from the panel.

**Mixture model.** Let $\mathbf{x}_i$ be the panel counts of cell $i$, with total $n_i$. The model is

$$\mathbf{x}_i \sim \operatorname{Multinomial}\!\left(n_i,\ \textstyle\sum_k w_{ik}\,\mathbf{p}_k\right),$$

where $\mathbf{p}_k$ is the fixed profile of subclass $k$ and $w_{ik}$ are per-cell mixing weights with $\sum_k w_{ik} = 1$.
- The profile $\mathbf{p}_k$ is the mean per-cell panel proportion of the subclass's reference cells, plus a pseudocount of 10⁻⁶, renormalised to sum to 1.
- The weights were estimated by expectation–maximisation (EM), starting from uniform weights and using the update

$$w_{ik} \leftarrow \frac{w_{ik}}{n_i} \sum_g \frac{x_{ig}\, p_{kg}}{\sum_l w_{il}\, p_{lg}}.$$

- EM stopped when no weight changed by more than 10⁻⁴, or after 300 iterations.

**Reference refinement.** References were refined twice. In each round:
1. Mixture weights were fitted for all reference cells.
2. Only the half of each subclass's reference cells with the highest weight on its own lineage was kept.
3. The profiles were recomputed from these cells.

This keeps the minority of contaminated cells from leaking into the profiles.

**Flagging.** With the final profiles, the foreign score of cell $i$ was $f_i = 1 - \sum_{k:\ \mathrm{lineage}(k) = \mathrm{lineage}(i)} w_{ik}$. This is the fraction of the cell's panel UMIs attributed to other lineages: about 0 for clean cells and about 0.5 for a doublet of two cells with similar panel depth. Cells with $f_i \ge 0.3$ were flagged and removed.

## Removal of embedding outliers

Among the cells not flagged as mixed, cells isolated from the rest of their subclass in the integrated embedding were identified with the local outlier factor (LOF; ref. 21; scikit-learn, via `sct.tl.detect_outliers`). LOF was fitted separately for each curated subclass, on the 30 Harmony-corrected components from the re-embedding in step 09, with Euclidean distance.
- Number of neighbours: 4% of the subclass's cells, bounded to the range 100–4,000 and at most *n* − 1.
- Cells with a LOF score above 1.5 were classified as outliers (scikit-learn `contamination="auto"`).
- Subclasses with fewer than 100 cells were not assessed.

Cells flagged as mixed or as outliers were removed on disk (`10_clean.py`), and UMAPs of the removed and retained cells were saved.

## Final embedding and release data

The retained cells were re-embedded as in step 09 (6,000 HVFs, PCA, regression, Harmony, kNN graph and UMAP; `11_umap.py`), and the result was written as the release H5AD file. It contains:
- raw counts (`raw.X`) and log-normalised expression (`X`)
- the curated class, subclass and subtype labels
- QC metrics, signature scores and doublet scores
- the corrected components and the UMAP coordinates

## Reproducibility

All stochastic steps used random seed 0: PCA, Harmony, neighbour-index construction, Leiden, UMAP, doublet simulation, the *k*-means steps of signature scoring, and reference sampling for the mixing model.

Two aspects can still vary between runs:
- With multithreading, Harmony can differ between runs at the order of 10⁻⁵ in the corrected components. UMAP coordinates can amplify these differences. Single-threaded Harmony (`n_jobs=1`) gives bit-identical results.
- Outputs depend on the software versions in Table 2. All steps were run in one conda environment (`py312`), orchestrated by Snakemake (ref. 22).

**Table 2 | Software and versions** (reference environment used to validate pipeline version 1.0.1).

| Software | Version | Use |
| --- | --- | --- |
| Python | 3.12.11 | runtime |
| SCTools | 0.5.0 | on-disk H5AD operations, HVF selection, QC thresholds, mixing and outlier detection |
| pegasuspy | 1.10.2 | QC metrics, PCA, regression, neighbours, Leiden, UMAP wrapper, Scrublet, signature scores |
| pegasusio | 0.9.1 | data input/output |
| scanpy | 1.11.0 | normalisation, HVF selection, plotting |
| anndata | 0.10.9 | H5AD data structures |
| harmony-pytorch | 0.1.8 | Harmony integration |
| PyTorch | 2.6.0 | Harmony back end |
| hnswlib | 0.8.0 | approximate nearest neighbours |
| leidenalg | 0.10.2 | Leiden clustering |
| igraph | 0.11.8 | graph back end |
| umap-learn | 0.5.7 | UMAP |
| pynndescent | 0.5.13 | umap-learn dependency |
| scikit-learn | 1.5.2 | PCA, regression, exact neighbours, *k*-means, LOF |
| statsmodels | 0.14.4 | MAD in HVF selection, Benjamini–Hochberg correction |
| NumPy | 1.26.4 | numerical computing |
| SciPy | 1.15.2 | sparse matrices, MAD in QC thresholds |
| pandas | 2.2.3 | tables |
| h5py | 3.13.0 | HDF5 access |
| matplotlib | 3.10.1 | figures |
| seaborn | 0.13.2 | QC figures |

**Table 3 | Curation software and versions** (recorded in the notebook logs of `milestone-3.3-taxonomy-integration`; CAVATICA Data Studio, up to 2 TB of memory).

| Software | Version | Use |
| --- | --- | --- |
| pyMN (pymetaneighbor) | 0.1.0, commit `882b545` of https://github.com/gillislab/pyMN | variable genes, MetaNeighbor AUROCs, top hits |
| scanpy | 1.11.5 | data input/output |
| anndata | 0.12.7–0.12.18 | H5AD data structures |
| NumPy | 1.26.4 | numerical computing |
| pandas | 2.1.1 | tables |
| SciPy | 1.17.0–1.17.1 | sparse matrices |
| R, with MetaNeighbor, dplyr, tidyr and ComplexHeatmap | not recorded | consistency filters, curation tables, figures |

**Reference data.**
- Ensembl release 104 gene annotation for GRCh38 (`annotation/all_genes_EnsemblInfo_104.csv`; ref. 23)
- Human MitoCarta3.0 (`annotation/Human.MitoCarta3.0.csv`; ref. 1)
- GENESIS colour palette (`annotation/260122_GENESIS_color_palette.csv`)

## Code availability

The pipeline code is in this directory (GENESIS unified processing pipeline v1.0.1, forked from GENESIS_iClust v1.0.1). SCTools is available at https://github.com/hoondy/SCTools (v0.5.0). The curation code is in [`milestone-3.3-taxonomy-integration`](../milestone-3.3-taxonomy-integration).

## References

1. Rath, S. et al. MitoCarta3.0: an updated mitochondrial proteome now with sub-organelle localization and pathway annotations. *Nucleic Acids Res.* **49**, D1541–D1547 (2021).
2. Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol.* **19**, 15 (2018).
3. Li, B. et al. Cumulus provides cloud-based data analysis for large-scale single-cell and single-nucleus RNA-seq. *Nat. Methods* **17**, 793–798 (2020).
4. Tirosh, I. et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science* **352**, 189–196 (2016).
5. Jerby-Arnon, L. et al. A cancer cell program promotes T cell exclusion and resistance to checkpoint blockade. *Cell* **175**, 984–997 (2018).
6. Zheng, G. X. Y. et al. Massively parallel digital transcriptional profiling of single cells. *Nat. Commun.* **8**, 14049 (2017).
7. Pedregosa, F. et al. Scikit-learn: machine learning in Python. *J. Mach. Learn. Res.* **12**, 2825–2830 (2011).
8. Korsunsky, I. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat. Methods* **16**, 1289–1296 (2019).
9. Malkov, Y. A. & Yashunin, D. A. Efficient and robust approximate nearest neighbor search using hierarchical navigable small world graphs. *IEEE Trans. Pattern Anal. Mach. Intell.* **42**, 824–836 (2020).
10. Haghverdi, L., Buettner, F. & Theis, F. J. Diffusion maps for high-dimensional single-cell analysis of differentiation data. *Bioinformatics* **31**, 2989–2998 (2015).
11. Coifman, R. R. & Lafon, S. Diffusion maps. *Appl. Comput. Harmon. Anal.* **21**, 5–30 (2006).
12. Traag, V. A., Waltman, L. & van Eck, N. J. From Louvain to Leiden: guaranteeing well-connected communities. *Sci. Rep.* **9**, 5233 (2019).
13. McInnes, L., Healy, J. & Melville, J. UMAP: uniform manifold approximation and projection for dimension reduction. Preprint at https://arxiv.org/abs/1802.03426 (2018).
14. Wolock, S. L., Lopez, R. & Klein, A. M. Scrublet: computational identification of cell doublets in single-cell transcriptomic data. *Cell Syst.* **8**, 281–291 (2019).
15. Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* **57**, 289–300 (1995).
16. Crow, M., Paul, A., Ballouz, S., Huang, Z. J. & Gillis, J. Characterizing the replicability of cell types defined by single cell RNA-sequencing data using MetaNeighbor. *Nat. Commun.* **9**, 884 (2018).
17. Fischer, S., Crow, M., Harris, B. D. & Gillis, J. Scaling up reproducible research for single-cell transcriptomics using MetaNeighbor. *Nat. Protoc.* **16**, 4031–4067 (2021).
18. Gabitto, M. I. et al. Integrated multimodal cell atlas of Alzheimer's disease. *Nat. Neurosci.* **27**, 2366–2383 (2024).
19. Ma, S. et al. Molecular and cellular evolution of the primate dorsolateral prefrontal cortex. *Science* **377**, eabo7257 (2022).
20. Siletti, K. et al. Transcriptomic diversity of cell types across the adult human brain. *Science* **382**, eadd7046 (2023).
21. Breunig, M. M., Kriegel, H.-P., Ng, R. T. & Sander, J. LOF: identifying density-based local outliers. In *Proc. 2000 ACM SIGMOD International Conference on Management of Data* 93–104 (2000).
22. Mölder, F. et al. Sustainable data analysis with Snakemake. *F1000Research* **10**, 33 (2021).
23. Howe, K. L. et al. Ensembl 2021. *Nucleic Acids Res.* **49**, D884–D891 (2021).
