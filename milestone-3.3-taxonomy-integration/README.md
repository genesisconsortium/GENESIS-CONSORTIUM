# GENESIS cortical subclass annotation

This repository contains interactive Python notebooks and R scripts used to build, inspect, and validate harmonized cortical cell-subclass annotations with MetaNeighbor. The files were designed to run in the CAVATICA Data Studio environment, with required packages installed in each notebook or script and imports kept to a minimum. They do not form a linear pipeline; analyses were rerun and reviewed iteratively using cross-dataset top hits, AUROC scores, summary tables, and visualizations.

## Python notebooks

- `merge_human_cortex_datasets.ipynb` assembles a downsampled prefrontal cortex reference from SEA-AD, Ma et al., Siletti et al., PsychAD/GEN_A1, and AMP PD/GEN_A3, retaining genes shared across datasets.
- `Metaneighbor_largemem.ipynb` compares GEN_A1, GEN_A2, or cortical GEN_A3 subtypes with the five-dataset reference using the Python implementation of MetaNeighbor and writes AUROC and top-hit results. These large expression-matrix comparisons can require substantial RAM; the largest analyses were run on a machine with 2 TB of memory.
- `Metaneighbor_GEN_A1-3.cortex.subclass.ipynb` compares target datasets with the harmonized GEN_A1–A3 cortical subclass reference and generates all-vs-all, one-vs-best, and thresholded top-hit results.
- `Metaneighbor_GEN_A1-2.cortex.subclass.ipynb` compares region-specific GEN_A3 cortical datasets and GEN_A6 with the harmonized GEN_A1–A2 subclass reference after GEN_A3 was split by region.
- `GEN_A_all.cortex.subclass.ipynb` combines the harmonized GEN_A1–A3 annotations with annotated target datasets and reruns MetaNeighbor to assess subclass agreement across the full cortical collection.

## R scripts

- `GEN_subclass_python_metaneighbor_plots.R` interactively summarizes and visualizes MetaNeighbor results, applies the cross-dataset consistency filters, and writes passing, failing, and assigned-subclass tables.
- `combine_collapsed.R` combines the GEN_A1–A3 top-hit profiles used to review consistent cross-reference matches and prepare the harmonized subclass annotation table.
