#!/bin/bash
#SBATCH --job-name=fnxfig
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
#SBATCH --output=fornax_figure_%j.out
#SBATCH --error=fornax_figure_%j.err

source /etc/profile.d/z00_lmod.sh
module purge
module load python/3.9.1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg
export MPLCONFIGDIR="$PWD/.mplconfig"

python3 -c 'import h5py, numpy, scipy, matplotlib, astropy; print("figure environment ready")'
python3 make_evolution_montage.py /travail/xhuang/Fornax2073 \
  --elinfo /travail/xhuang/Fornax2073/elinfo_Fornax2073.csv \
  --lookback-step-gyr 0.5 \
  --outdir output
