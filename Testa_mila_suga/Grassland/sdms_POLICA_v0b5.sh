#!/bin/bash
#SBATCH --job-name=POLICA_v0b5		 # Job name
#SBATCH --partition=power			 # Partition name
#SBATCH --ntasks=1				 # Number of tasks
#SBATCH --cpus-per-task=1
#SBATCH --mem=35G				 # Kopējā atmiņa
#SBATCH --time=72:00:00			 # Time limit, hrs:min:sec
#SBATCH --output=TestingScripts/JekaterinaButkevica/Papilionoidea/00GrasslandIndex/sdms_POLICA_v0b5.out			 # Standard output and error log



set -euo pipefail

# Absolute paths (avoid ../ relative surprises)
PROJECT_ROOT="/home/hiqbiodiv"
IMG="${PROJECT_ROOT}/hiqbiodiv-container_20260301.sif"
WORKDIR="${PROJECT_ROOT}"
SCRIPT="${WORKDIR}/TestingScripts/JekaterinaButkevica/sdms_POLICA_v0b5.R"

echo "Starting job"
echo "Date = $(date)"
echo "Node hostname = $(hostname -s)"
echo "Submit dir  = $(pwd)"

module load singularity

echo "Host check:"
ls -l "$SCRIPT" || { echo "Script not found on host: $SCRIPT"; exit 1; }

# Bind the project into the container and set the PWD inside the container
# Bind to the same path to keep file references identical.
echo "Running in container..."
srun singularity exec \
  -B "${PROJECT_ROOT}:${PROJECT_ROOT}" \
  --pwd "${WORKDIR}" \
  "${IMG}" \
  Rscript --vanilla "${SCRIPT}"

echo "All done on $(date)"