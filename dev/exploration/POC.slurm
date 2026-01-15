#!/bin/bash
#SBATCH --account=r01304
#SBATCH --mail-user=zzahid@iu.edu
#SBATCH --partition=general
#SBATCH --job-name=cm_poc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=128G
#SBATCH --time=5:00:00
#SBATCH --array=1-20
#SBATCH --output=/dev/null

export HOME="/N/u/zzahid/Quartz"

LOG_DIR="$HOME/Documents/Github/hrf-z/dev/exploration/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/cm_poc_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

echo "[$(date)] Starting array task $SLURM_ARRAY_TASK_ID" | tee -a "${LOG_FILE}"

module purge
module load StdEnv
module load gnu/12.2.0
module load r/4.4.1
module load rstudio/2023.09.1
module load java/15.0.2
module load mvapich/2.3.5
module load zlib/1.2.13
module load udunits/2.2.28
module load hdf5/1.12
module load netcdf/4.7
module load proj/8.0.1
module load gdal/3.6.4
module load gsl/2.6
module load sqlite/3.35.5
module load geos/3.8.2

Rscript "$HOME/Documents/Github/hrf-z/dev/exploration/POC_slurm.R" "$SLURM_ARRAY_TASK_ID" >> "${LOG_FILE}" 2>&1
rc=$?

if [ $rc -ne 0 ]; then
    echo "[$(date)] Task $SLURM_ARRAY_TASK_ID FAILED with exit code $rc" | tee -a "${LOG_FILE}"
    exit $rc
fi

echo "[$(date)] Task $SLURM_ARRAY_TASK_ID completed successfully" | tee -a "${LOG_FILE}"