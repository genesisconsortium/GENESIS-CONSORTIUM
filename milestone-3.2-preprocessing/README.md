# GENESIS whole-genome sequencing pipeline

This repository contains the scripts used to process whole-genome sequencing (WGS) data for the GENESIS Freeze v1 release. The workflow starts with paired-end FASTQ files, aligns the reads to the GRCh38 human reference genome, produces analysis-ready BAMs and per-sample genomic VCFs (gVCFs), performs cohort-level joint genotyping, applies variant quality score recalibration, and generates the files and quality-control metrics used for the final release.

The pipeline was run on an LSF high-performance computing cluster. Most stages were submitted as job arrays, with one task per sample, chromosome, or genomic interval. The scripts are numbered in execution order.

> **Important:** Paths, LSF project/queue names, array sizes, reference files, and cohort-specific sample manifests must be configured for the target environment before the pipeline is submitted.

## Workflow overview

| Step | Script | Main operation | Principal output |
|---:|---|---|---|
| 0 | `00_fastqc.sh` | Raw-read quality assessment | FastQC reports |
| 1 | `01_bwa.sh` | Alignment to GRCh38 with BWA-MEM | Aligned SAM/BAM |
| 2 | `02_samtools.sh` | BAM conversion, sorting, and indexing | Coordinate-sorted BAM |
| 3 | `03_fixmate.sh` | Mate-information correction and duplicate processing | Duplicate-marked BAM |
| 4 | `04_bqsr.sh` | GATK base quality score recalibration | Recalibrated BAM |
| 5 | `05_haplotypecaller.sh` | Per-sample variant calling in GVCF mode | Per-sample/per-chromosome gVCFs |
| 6 | `06_genomicsdbimport.sh` | Import gVCFs into interval-specific GenomicsDB workspaces | GenomicsDB workspaces |
| 7 | `07_genotypegvcfs.sh` | Joint genotyping of each genomic interval | Interval-level VCF/BCF |
| 8 | `08_concat.sh` | Trim padded intervals and concatenate results | Chromosome-level callsets |
| 9 | `09_vqsr.sh` | Separate SNP and indel VQSR | Recalibrated callsets |
| 10 | `10_finalfiles.sh` | Normalize, index, checksum, and assemble release files | Final joint-called VCF/BCF |
| 11 | `11_sample_QC.sh` | Sample-level sequencing and genotype QC | Sample QC tables |
| 12 | `12_variant_QC.sh` | Cohort- and variant-level QC | Variant QC tables and plots |

The repository may group the first five stages under **STOCHEUSIS** (alignment and per-sample processing), the joint-calling and cohort-QC stages under **SYNTHESIS**, and downstream functional annotation under **SEMANSIS**.

## Input data and required files

The workflow requires:

- paired-end FASTQ files for every sample;
- a tab-delimited sample manifest containing at least the sample identifier and the paths to read 1 and read 2;
- the GRCh38 reference FASTA and its BWA, SAMtools, and Picard/GATK index files;
- known-sites VCFs for BQSR;
- truth, training, and annotation resources for SNP and indel VQSR;
- chromosome and interval lists compatible with the reference contig names;
- sufficient temporary and permanent storage for BAMs, gVCFs, GenomicsDB workspaces, and joint-called files.

All reference resources must use the same genome build and contig naming convention. Mixing `chr1`-style and `1`-style contig names will cause interval, normalization, and concatenation failures.

An example manifest is:

```text
sample_id\tfastq_r1\tfastq_r2
SAMPLE001\t/path/to/SAMPLE001_R1.fastq.gz\t/path/to/SAMPLE001_R2.fastq.gz
SAMPLE002\t/path/to/SAMPLE002_R1.fastq.gz\t/path/to/SAMPLE002_R2.fastq.gz
```

## Software

The Freeze v1 implementation used the following principal tools:

- FastQC for read-level QC;
- BWA-MEM for alignment;
- SAMtools for BAM manipulation and indexing;
- Picard for duplicate marking and selected BAM/VCF metrics;
- GATK 4.6.1.0 for BQSR, HaplotypeCaller, GenomicsDBImport, GenotypeGVCFs, and VQSR;
- bcftools/htslib 1.21 for querying, trimming, concatenation, normalization, indexing, and final file manipulation;
- VerifyBamID/VerifyBamID2 for contamination assessment;
- PLINK/PLINK2 and KING for genotype-level sample QC, relatedness, ancestry, and sex checks.

Exact versions should be recorded in the run logs and release metadata. On the Mount Sinai LSF environment, software was loaded through environment modules within each job script.

## Directory structure

The analysis was organized approximately as follows:

```text
PROJECT/
├── 1_bwamem/       # initial alignments
├── 2_fixmate/      # mate-corrected/intermediate BAMs
├── 3_markdup/      # duplicate-marked BAMs
├── 4_BQSR/         # recalibrated BAMs
├── 5_GVCF/         # per-sample gVCFs, split by chromosome
├── 6_metrics/      # alignment and sample metrics
├── 7_joint/        # GenomicsDB and interval-genotyped files
├── 8_concat/       # chromosome-level concatenated callsets
├── 9_VQSR/         # VQSR models and recalibrated callsets
├── 10_final/       # final release files
├── logs/           # LSF stdout and stderr
├── scripts/        # pipeline scripts and configuration
├── tmp/            # temporary files
└── samples.txt     # ordered sample manifest
```

Each script creates its output and log directories when they do not already exist.

## Detailed pipeline

### 0. Raw-read QC (`00_fastqc.sh`)

FastQC was run on both FASTQ files from every sample before alignment. This stage was used to identify problems such as low per-base sequence quality, adapter contamination, abnormal GC content, overrepresented sequences, sequence duplication, and truncated or corrupted input files.

FastQC is diagnostic: reads were not discarded solely because a module was marked as a warning or failure. Reports were reviewed across the cohort to distinguish recurrent sequencing-platform patterns from sample-specific failures. The paired FASTQ files and sample identifiers were also checked against the manifest before computationally expensive alignment jobs were launched.

Typical submission:

```bash
bsub < scripts/00_fastqc.sh
```

For an array implementation, the LSF array index selects one row of the sample manifest and FastQC is run on both mates.

### 1. Alignment to GRCh38 (`01_bwa.sh`)

Reads were aligned to the harmonized GRCh38 reference with BWA-MEM. Each alignment included a read-group string containing, at minimum, a unique read-group identifier, sample identifier, library identifier, platform unit, and sequencing platform. Correct read groups are essential because GATK uses the `SM` field to associate reads with a biological sample.

Conceptually, this stage performs:

```bash
bwa mem \
  -R '@RG\tID:<RG_ID>\tSM:<SAMPLE>\tLB:<LIBRARY>\tPL:ILLUMINA\tPU:<PLATFORM_UNIT>' \
  reference.fa sample_R1.fastq.gz sample_R2.fastq.gz
```

The actual script supplies the project-specific thread count, reference location, read-group values, and output paths. LSF stdout and stderr are retained for every sample so alignment failures can be traced without relying only on the existence of an output BAM.

### 2. BAM conversion, sorting, and indexing (`02_samtools.sh`)

The BWA output was converted to BAM and coordinate-sorted. Sorting places alignments in reference-coordinate order, which is required by downstream Picard and GATK stages. The sorted BAM was indexed and subjected to a quick integrity check.

Typical operations include:

```bash
samtools view -b <alignment.sam> |
  samtools sort -o <sample.sorted.bam> -
samtools index <sample.sorted.bam>
samtools quickcheck -v <sample.sorted.bam>
```

Where the alignment script streams BWA output directly into SAMtools, the intermediate SAM is omitted. This reduces storage and I/O without changing the resulting coordinate-sorted BAM.

### 3. Mate correction and duplicate marking (`03_fixmate.sh`)

Mate fields were corrected before duplicate processing, and PCR/optical duplicates were marked. Duplicate marking retains the reads but flags redundant fragments so that GATK does not treat them as independent observations during variant calling.

Depending on the exact implementation, mate correction and duplicate marking may be performed with SAMtools (`collate`/name sort, `fixmate`, coordinate sort, and `markdup`) or with Picard. Freeze v1 release documentation records Picard duplicate marking as the authoritative duplicate-processing step. The duplicate-marking metrics file was retained for sample QC.

The final product of this stage is a coordinate-sorted, indexed, duplicate-marked BAM. BAM integrity, header/sample identity, sort order, and index availability were checked before BQSR.

### 4. Base quality score recalibration (`04_bqsr.sh`)

GATK BaseRecalibrator was run against GRCh38 and approved known-variant resources to model systematic errors in the sequencer-assigned base quality scores. GATK ApplyBQSR then applied the recalibration model to the duplicate-marked BAM.

Conceptually:

```bash
gatk BaseRecalibrator \
  -R reference.fa \
  -I sample.markdup.bam \
  --known-sites known_sites_1.vcf.gz \
  --known-sites known_sites_2.vcf.gz \
  -O sample.recal.table

gatk ApplyBQSR \
  -R reference.fa \
  -I sample.markdup.bam \
  --bqsr-recal-file sample.recal.table \
  -O sample.bqsr.bam
```

The recalibrated BAM was indexed and used as the definitive read-level input for variant calling and final sample QC. BQSR tables and logs were retained as provenance.

### 5. Per-sample gVCF generation (`05_haplotypecaller.sh`)

GATK HaplotypeCaller 4.6.1.0 was run independently for each sample in reference-confidence mode:

```bash
gatk HaplotypeCaller \
  -R reference.fa \
  -I sample.bqsr.bam \
  -O sample.chrN.g.vcf.gz \
  -L chrN \
  --emit-ref-confidence GVCF
```

Calling was partitioned by chromosome (`1`-`22`, `X`, `Y`, and other reference contigs as configured) to distribute the workload and simplify reruns. Each output was bgzip-compressed, tabix-indexed, and accompanied by an MD5 checksum. A complete sample therefore produced one gVCF, one `.tbi`, and one checksum per configured chromosome partition.

Before joint calling, the pipeline verified that every expected sample/partition combination existed, was nonempty, passed a basic format check, had a readable index, and used sample and contig names consistent with the manifest and reference.

### 6. GenomicsDB import (`06_genomicsdbimport.sh`)

Joint calling was divided into approximately 1-Mb nonoverlapping **core intervals**. Each core interval was expanded by approximately 10 kb on both sides for import and genotyping. Padding provides local haplotype context at interval boundaries; it is removed before final concatenation so each genomic position appears only once.

For each padded interval, the per-sample gVCFs were imported into a GATK GenomicsDB workspace. At the scale of Freeze v1, imports were performed in batches of approximately 250 samples to manage file-handle, memory, and runtime limits. Sample-name maps were generated from the ordered manifest and contained one unique sample-to-gVCF mapping per input.

Conceptually:

```bash
gatk GenomicsDBImport \
  --genomicsdb-workspace-path <workspace> \
  --sample-name-map <sample_map.tsv> \
  --intervals <padded_interval> \
  --reader-threads <N>
```

GenomicsDB workspaces are directories rather than single files. A successful job was determined from the exit status and log, not merely from the presence of the workspace directory. Failed or partial workspaces were rebuilt for the exact affected interval/batch rather than reused.

### 7. Joint genotyping (`07_genotypegvcfs.sh`)

GATK GenotypeGVCFs was run on each interval-specific GenomicsDB workspace to convert reference-confidence records into jointly genotyped cohort variants:

```bash
gatk GenotypeGVCFs \
  -R reference.fa \
  -V gendb://<workspace> \
  -L <padded_interval> \
  -O <interval>.vcf.gz
```

Joint genotyping was performed independently for all intervals. This approach makes the cohort-wide computation restartable: a failed interval can be rerun without regenerating an entire chromosome. Interval outputs were indexed and validated before concatenation.

### 8. Core-interval trimming and chromosome concatenation (`08_concat.sh`)

Each jointly genotyped interval was restricted to its original, nonoverlapping core coordinates. This removes duplicate calls introduced by the 10-kb padding. Core interval files were ordered numerically by genomic coordinate and concatenated with bcftools to produce one callset per chromosome.

Conceptually:

```bash
bcftools view -r <core_interval> <padded_interval.vcf.gz> -Oz -o <core.vcf.gz>
bcftools concat -f <ordered_core_file_list> -Ob -o Joint_chrN.bcf
bcftools index Joint_chrN.bcf
```

The ordered input list was checked for missing intervals, accidental duplication, inconsistent sample order, overlapping cores, and gaps relative to the expected interval definition. Concatenation should never be performed using lexicographic filename order unless filenames are zero-padded; otherwise interval 100 can incorrectly precede interval 20.

### 9. Variant quality score recalibration (`09_vqsr.sh`)

GATK VariantRecalibrator and ApplyVQSR were run separately for SNPs and indels using established truth, training, and annotation resources compatible with GRCh38. Separate modeling is required because SNPs and indels have different error profiles.

The general sequence was:

1. build the SNP recalibration model;
2. apply the selected SNP truth-sensitivity tranche;
3. build the indel recalibration model on the SNP-recalibrated callset;
4. apply the selected indel tranche;
5. retain the VQSR model files, tranche files, plots, logs, and resource definitions.

Variants outside the accepted tranche remain represented with their FILTER status unless a release-specific downstream step explicitly selects `PASS` records. The chosen tranche thresholds and resource priors are scientific parameters and should be recorded in the script/configuration and release metadata.

For smaller chromosome partitions or variant classes that do not contain enough observations to fit a stable VQSR model, the pipeline should use a documented cohort-appropriate strategy rather than silently treating an unsuccessful model as passing QC.

### 10. Final callset assembly (`10_finalfiles.sh`)

The finalization stage standardized and validated the release files. Operations included, as applicable:

- normalization and left-alignment against GRCh38;
- checking REF alleles against the reference FASTA;
- ensuring consistent headers, contig declarations, sample names, and sample order;
- indexing every final VCF/BCF;
- concatenating chromosome-level outputs into the intended release layout;
- generating MD5 checksums;
- recording sample manifests, software versions, commands, and run provenance;
- running Picard and bcftools summary metrics.

Typical validation commands include:

```bash
bcftools index -s final.vcf.gz
bcftools query -l final.vcf.gz
bcftools stats final.vcf.gz > final.vcf.stats
md5sum final.vcf.gz final.vcf.gz.tbi > final.md5
```

The release is considered complete only when the expected number and order of samples, contig coverage, indexes, checksums, FILTER fields, and summary counts have been verified.

## Quality control

QC was performed throughout the pipeline rather than only after joint calling.

### Alignment and sample QC (`11_sample_QC.sh`)

Read- and BAM-level metrics included:

- total reads and alignment rate;
- properly paired read fraction;
- duplicate fraction;
- insert-size distribution;
- mean and median depth;
- breadth of coverage at selected thresholds;
- base and mapping quality summaries;
- BAM integrity and read-group/sample consistency;
- contamination estimates from VerifyBamID/VerifyBamID2.

Historically, samples with `FREEMIX > 0.03` were flagged for contamination review. This is a review threshold, not a substitute for examining the complete evidence and cohort context.

Genotype-based sample QC included:

- per-sample missingness/call rate;
- mean depth and genotype quality;
- heterozygosity and heterozygous/homozygous-alternate ratios;
- transition/transversion ratio;
- genetic-sex inference and reported-sex concordance;
- chrX and chrY call-rate and coverage checks;
- identity and duplicate detection;
- pairwise relatedness with KING/PLINK;
- ancestry inference and PCA;
- sample outlier review across sequencing batches and data sources.

For Freeze v1, the release-level QC documentation also tracked samples excluded from the final analytic set and preserved intentional duplicates when required for concordance assessment.

### Variant and cohort QC (`12_variant_QC.sh`)

Variant-level QC included:

- total SNP and indel counts;
- PASS and filtered counts by chromosome and variant type;
- multiallelic-site counts;
- missingness/call rate;
- allele count and allele-frequency distributions;
- singleton counts;
- depth and genotype-quality distributions;
- Ti/Tv ratios;
- heterozygous/homozygous-alternate ratios;
- Hardy-Weinberg equilibrium in an appropriate unrelated population subset;
- VQSR tranche and annotation distributions;
- concordance or replicate checks where duplicate data were available;
- sex-chromosome-specific QC, including PAR/non-PAR handling and chrY evaluation in genetically inferred males.

QC summaries were reviewed by chromosome, ancestry, sequencing batch, source cohort, and case/control or phenotype group where appropriate. Stratification helps distinguish technical artifacts from expected biological differences.

## Running the workflow

1. Edit the central configuration or script variables to define the project directory, sample manifest, reference FASTA, resource VCFs, interval files, LSF project/queue, memory, wall time, and software modules.
2. Validate the manifest and confirm that all FASTQs are readable.
3. Submit each numbered stage only after the previous stage has completed successfully.
4. Review exit codes, stderr, expected file counts, indexes, and validation metrics before advancing.
5. Rerun only failed array elements or genomic intervals when possible.

A simple submission pattern is:

```bash
bsub < scripts/00_fastqc.sh
bsub < scripts/01_bwa.sh
bsub < scripts/02_samtools.sh
bsub < scripts/03_fixmate.sh
bsub < scripts/04_bqsr.sh
bsub < scripts/05_haplotypecaller.sh
bsub < scripts/06_genomicsdbimport.sh
bsub < scripts/07_genotypegvcfs.sh
bsub < scripts/08_concat.sh
bsub < scripts/09_vqsr.sh
bsub < scripts/10_finalfiles.sh
bsub < scripts/11_sample_QC.sh
bsub < scripts/12_variant_QC.sh
```

In production, job dependencies or explicit completion checks should prevent a downstream stage from starting when upstream array elements are missing or failed. Do not infer successful completion solely from the scheduler reporting that a job has left the queue.

## Reproducibility and provenance

For every release, retain:

- the exact sample manifest and its checksum;
- the reference FASTA and resource versions/checksums;
- chromosome and padded/core interval definitions;
- complete job scripts and configuration;
- software/module versions;
- LSF stdout/stderr logs;
- VQSR recalibration, tranche, and plot files;
- final indexes and MD5 checksums;
- sample exclusion and inclusion tables;
- QC reports and release notes.

Changes in samples, reference resources, intervals, tool versions, or calling parameters define a new run and should be documented rather than overwriting the provenance of a frozen release.

## Known implementation considerations

- **Interval padding:** padding is used for computational context only; padded records must be trimmed back to nonoverlapping core intervals before concatenation.
- **Sample order:** all interval and chromosome files must contain the same samples in the same order.
- **Duplicate sample identifiers:** sample names must be unique within a joint-calling run. Replicates should receive unambiguous identifiers and be reconciled deliberately during QC.
- **Sex chromosomes:** chrX PAR/non-PAR and chrY require ploidy-aware interpretation; autosomal QC rules cannot always be applied unchanged.
- **Alternative contigs:** `other` contigs must be defined explicitly and handled consistently across HaplotypeCaller, joint calling, and final concatenation.
- **GenomicsDB workspaces:** a partially created workspace is not a valid checkpoint. Inspect logs and rebuild failed intervals.
- **Storage:** intermediate BAMs, chromosome-partitioned gVCFs, and GenomicsDB directories are storage-intensive. Check available capacity before each major stage.
- **Checksums:** an MD5 for the gVCF does not validate its `.tbi`; both files must exist and be readable, and release artifacts should be checksummed explicitly.

## Citation and acknowledgement

When using this workflow or its callsets, cite the appropriate GENESIS data release, reference genome, software packages, and contributing cohorts. Add project-specific acknowledgement and data-access language here before public release.

## Contact

For questions about the pipeline or Freeze v1 release, open a GitHub issue or contact the GENESIS analysis team.
