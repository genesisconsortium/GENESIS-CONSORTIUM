# GENESIS-CONSORTIUM
# GENESIS-CONSORTIUM 
**Gene Expression aNd Epigenetics Single cell Integrative Study** *Constructing a cell type-specific quantitative trait locus (QTL) atlas for the human brain.*

---

## 🧬 Program Overview
GENESIS (BrainCellQTL) addresses the critical need for cell-type-specific biological mechanisms in neuropsychiatric and neurodegenerative disorders (Alzheimer’s, Parkinson’s, Schizophrenia). By harmonizing over **10,000 single-cell libraries** from **3,000+ unique brain donors**, we are constructing a comprehensive xQTL atlas to:

1.  **Increase mechanistic understanding** of cellular dysfunction in these disorders.
2.  **Better prioritize significant genes and molecular pathways** for future hypothesis-driven mechanistic studies.
3.  **Provide a valuable resource** that can be applied in ongoing and future genome-wide association studies.
4.  **Provide preprocessed and harmonized single-nucleus brain omics data** for the research community.
5.  **Establish a mechanism for data sharing and harmonization** across consortia for future mega-analyses.

---

## 🧭 User Flow Navigation
*Select your path below to access relevant documentation and pipelines.*

* **📊 Processed Data User:** [AIM3 Overview](./genesis-aim3-overview) → [Release Notes](./genesis-docs) → [Output Specs](./milestone-3.2-preprocessing/single-cell/shared)
* **🧬 Pipeline User:** [Select Milestone](#-repository-map) → [Cloud (CAVATICA)](./) or [Local (HPC)](./) → [README](./)
* **💻 Method Developer:** [Local Setup](./) → [Docs](./genesis-docs) → [Assumptions & Limitations](./milestone-3.2-preprocessing/docs)

---

## 📂 Repository Map

### [3.1 — Establish Pipelines on CAVATICA](./milestone-3.1-pipelines)
*Infrastructure for cloud-native and local high-performance computing.*
* **Cloud:** CWL & Dockerized workflows.
* **Local:** Snakemake and Nextflow implementations for HPC environments.

### [3.2 — Preprocess Single-Cell & Genotype Data](./milestone-3.2-preprocessing)
*Standardized preprocessing for multi-modal datasets.*
* **Single-Cell:** scRNAseq, scATACseq, and Multiome (10x) pipelines.
* **Genotype:** WGS (Joint-calling) and TR-VNTR-SV pipelines.

### [3.3 — Taxonomy & Cross-Cohort Integration](./milestone-3.3-taxonomy-integration)
*Reference taxonomy generation and batch-effect correction.*
* **Reference Taxonomy:** Building the "Zombosome" (mitochondrial) and cell-type atlas.
* **Integration:** MetaNeighbor and MetaMarkers for multi-site harmonization.

### [3.4 — QTL & TR-QTL Analyses](./milestone-3.4-qtl-analysis)
*Mapping genetic variants to molecular phenotypes.*
* **QTL:** cis-QTL models and exploratory trans-QTL analysis.
* **TR-QTL:** Specialized mapping for STRs, VNTRs, and Structural Variants.

### [3.5 — Genomic Feature Imputation (GFI)](./milestone-3.5-gfi-analysis)
*Predictive modeling and downstream integration.*
* **Models:** Feature sources and imputation training for gene dysregulation.

---

## Consortium Information
* **Institutions:** Icahn School of Medicine at Mount Sinai, Sage Bionetworks, University of Toronto, Velsera.
* **Grant Number:** U24AG087563
* **Contact:** [genesis.consortium@gmail.com](mailto:genesis.consortium@gmail.com)

---
© 2026 GENESIS Consortium | Icahn School of Medicine at Mount Sinai
