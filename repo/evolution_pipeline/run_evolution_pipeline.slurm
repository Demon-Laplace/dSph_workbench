#!/bin/bash
#SBATCH --job-name=fornax_evolution
#SBATCH --clusters=astro_thin
#SBATCH --partition=def
#SBATCH --account=yang
#SBATCH --qos=astro_thin_def_long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=fornax_evolution_%j.out
#SBATCH --error=fornax_evolution_%j.err

set -euo pipefail
: "${CONFIG_PATH:?Set CONFIG_PATH to an evolution-pipeline JSON configuration}"
: "${PIPELINE_ROOT:=/travail/xhuang/dSph_workbench/evolution_pipeline}"

. /etc/profile.d/z00_lmod.sh
module purge
module load python/3.9.1

export FORNAX_REPO_ROOT="${PIPELINE_ROOT}/.."
export MPLBACKEND=Agg
export MPLCONFIGDIR="${PIPELINE_ROOT}/.mplconfig"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
mkdir -p "${MPLCONFIGDIR}"

run_dir=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["paths"]["run_dir"])' "${CONFIG_PATH}")
cd "${run_dir}"
python3 "${PIPELINE_ROOT}/fornax_evolution_pipeline.py" run --config "${CONFIG_PATH}"
