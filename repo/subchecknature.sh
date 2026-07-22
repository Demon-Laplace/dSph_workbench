#!/bin/sh
#SBATCH --job-name=checknature
#SBATCH --mail-user=None
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --clusters=astro_thin
#SBATCH --partition=def
#SBATCH --qos=astro_thin_def_long
#SBATCH --account=yang
#SBATCH --time=01:00:00

module purge
module load python/3.9.1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

python3 ../dSph_workbench/PlotFig.py --processes "${SLURM_CPUS_PER_TASK:-24}"
