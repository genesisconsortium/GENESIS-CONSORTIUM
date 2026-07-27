#!/usr/bin/env bash
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -n 1
#BSUB -J "QC_RESID[1-15]"
#BSUB -W 02:00
#BSUB -M 16000
#BSUB -R "rusage[mem=16000]"
#BSUB -o /sc/arion/projects/CommonMind/genesis/pseudobulk_dreamlet_rc2/logs/QC_RESID_%I_%J.out
#BSUB -e /sc/arion/projects/CommonMind/genesis/pseudobulk_dreamlet_rc2/logs/QC_RESID_%I_%J.err

module purge
ml R/4.4.1

BASE=/sc/arion/projects/CommonMind/genesis/pseudobulk_dreamlet_rc2
CONFIG=${BASE}/config/pseudobulk_inputs.csv
SCRIPT=${BASE}/scripts/02_qc_pearson_residual_files.R

ROW=${LSB_JOBINDEX}

OUTPUT_DIR=$(awk -F',' -v row="${ROW}" 'NR == row + 1 {print $8}' "${CONFIG}")
RESIDUAL_DIR=${OUTPUT_DIR}/residuals

echo "LSB_JOBINDEX=${LSB_JOBINDEX}"
echo "CONFIG=${CONFIG}"
echo "OUTPUT_DIR=${OUTPUT_DIR}"
echo "RESIDUAL_DIR=${RESIDUAL_DIR}"
echo "SCRIPT=${SCRIPT}"

if [ ! -d "${RESIDUAL_DIR}" ]; then
  echo "ERROR: residual directory does not exist: ${RESIDUAL_DIR}"
  exit 1
fi

Rscript "${SCRIPT}" "${RESIDUAL_DIR}"
