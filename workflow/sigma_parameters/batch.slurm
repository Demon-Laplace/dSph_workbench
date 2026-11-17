#!/bin/bash
#SBATCH --clusters=mesopsl1
#SBATCH --partition=def
#SBATCH --qos=mesopsl1_def_long
#SBATCH --account=yang
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=56G
#SBATCH --time=99:00:00
#SBATCH --job-name=Fornax_sigma_data
#SBATCH --chdir=/travail/xhuang/workspace/jobs/sigma-parameters
#SBATCH --output=/travail/xhuang/workspace/jobs/sigma-parameters/compute-%j.out
#SBATCH --error=/travail/xhuang/workspace/jobs/sigma-parameters/compute-%j.err
set -e
source /etc/profile.d/z00_lmod.sh
module purge
module load python/3.9.1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 PYTHONIOENCODING=utf-8
export TMPDIR=/travail/xhuang/workspace/jobs/sigma-parameters/tmp
mkdir -p "$TMPDIR"
trap 'code=$?; if [ "$code" -ne 0 ]; then python3 -B /travail/xhuang/workspace/jobs/sigma-parameters/write_failure.py "$code"; fi' EXIT
python3 -B -u /travail/xhuang/workspace/jobs/sigma-parameters/prepare.py
python3 -B -u /travail/xhuang/workspace/jobs/sigma-parameters/run.py
