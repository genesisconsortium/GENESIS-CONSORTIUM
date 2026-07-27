#!/usr/bin/env bash
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -n 1
#BSUB -J "SEXQC[1,2,7,12,13,14,15]"
#BSUB -W 24:00
#BSUB -M 400000
#BSUB -R "rusage[mem=400000]"
#BSUB -o /sc/arion/projects/CommonMind/genesis/pseudobulk_dreamlet_rc2/logs/SEXQC_%I_%J.out
#BSUB -e /sc/arion/projects/CommonMind/genesis/pseudobulk_dreamlet_rc2/logs/SEXQC_%I_%J.err

set -euo pipefail

module purge
ml freetype libtiff libjpeg webp libxml2 binutils/2.38 R/4.5.1

BASE=/sc/arion/projects/CommonMind/genesis/pseudobulk_dreamlet_rc2
CONFIG=${BASE}/config/pseudobulk_inputs.csv
SEX_QC_FILE=${BASE}/sex_marker_qc/sex_qc_samples.csv
SCRIPT=${BASE}/scripts/03_rerun_pearson_residuals_sex_qc_pass.R

if [[ -z "${LSB_JOBINDEX:-}" ]]; then
  echo "ERROR: LSB_JOBINDEX is not set. Submit this script with bsub." >&2
  exit 1
fi

if [[ ! -f "${CONFIG}" ]]; then
  echo "ERROR: Config file not found: ${CONFIG}" >&2
  exit 1
fi

if [[ ! -f "${SEX_QC_FILE}" ]]; then
  echo "ERROR: Sex-QC file not found: ${SEX_QC_FILE}" >&2
  exit 1
fi

if [[ ! -f "${SCRIPT}" ]]; then
  echo "ERROR: R script not found: ${SCRIPT}" >&2
  exit 1
fi

# The LSF index refers to a data row in the config (header is not counted).
# Use R's CSV parser rather than shell field splitting so quoted CSV values are
# handled correctly.
PROJECT_OUTPUT_DIR=$(
  Rscript -e '
    config <- read.csv(
      commandArgs(trailingOnly = TRUE)[1],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    idx <- as.integer(commandArgs(trailingOnly = TRUE)[2])
    if (is.na(idx) || idx < 1 || idx > nrow(config)) {
      stop("LSB_JOBINDEX out of range: ", idx)
    }
    if (!"output_dir" %in% colnames(config)) {
      stop("Config is missing output_dir")
    }
    cat(config$output_dir[idx])
  ' "${CONFIG}" "${LSB_JOBINDEX}"
)

case "${LSB_JOBINDEX}" in
  1)  FILE_REGION=GEN_A1 ;;
  2)  FILE_REGION=GEN_A2 ;;
  7)  FILE_REGION=GEN_A4_MEC ;;
  12) FILE_REGION=GEN_A5 ;;
  13) FILE_REGION=GEN_A6 ;;
  14) FILE_REGION=GEN_A13 ;;
  15) FILE_REGION=GEN_A16 ;;
  *)
    echo "ERROR: No file_region configured for index ${LSB_JOBINDEX}" >&2
    exit 1
    ;;
esac

PROCESSASSAYS_DIR=${PROJECT_OUTPUT_DIR}/processassays
RERUN_OUTPUT_DIR=${PROJECT_OUTPUT_DIR}/sex_qc_pass_rerun
echo "FILE_REGION=${FILE_REGION}"

if [[ ! -d "${PROCESSASSAYS_DIR}" ]]; then
  echo "ERROR: processAssays directory not found: ${PROCESSASSAYS_DIR}" >&2
  exit 1
fi

mkdir -p "${RERUN_OUTPUT_DIR}" "${BASE}/logs"

echo "LSB_JOBINDEX=${LSB_JOBINDEX}"
echo "CONFIG=${CONFIG}"
echo "SEX_QC_FILE=${SEX_QC_FILE}"
echo "PROJECT_OUTPUT_DIR=${PROJECT_OUTPUT_DIR}"
echo "PROCESSASSAYS_DIR=${PROCESSASSAYS_DIR}"
echo "RERUN_OUTPUT_DIR=${RERUN_OUTPUT_DIR}"
echo "SCRIPT=${SCRIPT}"

Rscript "${SCRIPT}" \
  "${SEX_QC_FILE}" \
  "${PROCESSASSAYS_DIR}" \
  "${RERUN_OUTPUT_DIR}" \
  "${FILE_REGION}"

