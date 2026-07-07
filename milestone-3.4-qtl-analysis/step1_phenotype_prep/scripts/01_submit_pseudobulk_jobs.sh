#!/usr/bin/env bash
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -n 1
#BSUB -J "PB[12]"
#BSUB -W 24:00
#BSUB -M 915000
#BSUB -R "rusage[mem=915000]"
#BSUB -o /sc/arion/projects/CommonMind/genesis/pseudobulk_dreamlet_rc2/logs/PB_%I_%J.out
#BSUB -e /sc/arion/projects/CommonMind/genesis/pseudobulk_dreamlet_rc2/logs/PB_%I_%J.err

module purge
ml R/4.4.1
module load udunits
module load proj
module load geos
module load gdal

BASE=/sc/arion/projects/CommonMind/genesis/pseudobulk_dreamlet_rc2
CONFIG=${BASE}/config/pseudobulk_inputs.csv
SCRIPT=${BASE}/scripts/01_create_pseudobulk_residuals_from_h5ad.R

echo "LSB_JOBINDEX=${LSB_JOBINDEX}"
echo "CONFIG=${CONFIG}"
echo "SCRIPT=${SCRIPT}"

Rscript ${SCRIPT} ${CONFIG} ${LSB_JOBINDEX}
