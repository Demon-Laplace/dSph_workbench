#!/bin/sh
#SBATCH --job-name=IC_check
#SBATCH --mail-user=None    
#SBATCH --mail-type=ALL
#SBATCH --nodes=1                         
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24                                       
#SBATCH --clusters=astro_thin
#SBATCH --partition=hi
#SBATCH --qos=astro_thin_hi_normal
#SBATCH --account=yang_hi
#SBATCH --time=01:00:00

set -eu
cd "${SLURM_SUBMIT_DIR:-$PWD}"
isolation_root=$PWD
while [ ! -f "$isolation_root/python/check_density_profile.py" ]; do
    [ "$isolation_root" != / ] || { echo "Cannot locate isolation public Python" >&2; exit 1; }
    isolation_root=$(dirname "$isolation_root")
done

. /etc/profile.d/z00_lmod.sh
module purge
module load python/3.9.1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

exec python3 "$isolation_root/python/check_density_profile.py" --processes 24 --replace

