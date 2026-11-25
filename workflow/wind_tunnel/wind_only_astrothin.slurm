#!/bin/bash
#SBATCH --job-name=f8002_wind_nohost
#SBATCH --clusters=astro_thin
#SBATCH --partition=def
#SBATCH --qos=astro_thin_def_long
#SBATCH --account=yang
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=72000M
#SBATCH --time=99:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
source /etc/profile.d/z00_lmod.sh
module purge
module load intel/18.0.1 intelmpi/18.0.1 hdf5/1.8.21 python/3.9.1
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export I_MPI_FABRICS=shm I_MPI_HYDRA_BOOTSTRAP=slurm
export I_MPI_PIN=1 I_MPI_PIN_DOMAIN=core I_MPI_DEBUG=4
export PYTHONIOENCODING=utf-8 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
exec python3 tools/run_server_wind.py run --config run.json --mode "${1:-2}"
