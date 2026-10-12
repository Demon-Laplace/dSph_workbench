#!/bin/bash
#SBATCH --job-name=fnxfig2
#SBATCH --mail-user=None
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --clusters=astro_thin
#SBATCH --partition=def
#SBATCH --qos=astro_thin_def_long
#SBATCH --account=yang
#SBATCH --time=01:00:00
#SBATCH --output=fornax_figure2_%j.out
#SBATCH --error=fornax_figure2_%j.err

source /etc/profile.d/z00_lmod.sh
module purge
module load python/3.9.1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg
export MPLCONFIGDIR="${TMPDIR:-/tmp}/fornax_figure2_mplconfig"

model_dir=/travail/xhuang/Fornax2073_10mfd_RPF0.8
script=/travail/xhuang/dSph_workbench_paper/repo/paper_figures/fornax_observable_comparison/make_current_observables.py
output_dir=/travail/xhuang/Fornax2073_paper_figures/output/figure2

cd "$model_dir"
python3 "$script" "$model_dir" --output-dir "$output_dir"
