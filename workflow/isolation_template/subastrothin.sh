#!/bin/sh
#SBATCH --job-name=isolation 
#SBATCH --time=99:00:00
#SBATCH --mail-user=None
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --clusters=astro_thin
#SBATCH --partition=def
#SBATCH --qos=astro_thin_def_long
#SBATCH --account=yang



#export OMP_NUM_THREADS=1
#export I_MPI_COMPATIBILITY=4
#export LD_LIBRARY_PATH=$	LD_LIBRARY_PATH:/usr/lib64

set -eu
cd "${SLURM_SUBMIT_DIR:-$PWD}"
isolation_root=$PWD
while [ ! -f "$isolation_root/python/check_density_profile.py" ]; do
    [ "$isolation_root" != / ] || { echo "Cannot locate isolation public Python" >&2; exit 1; }
    isolation_root=$(dirname "$isolation_root")
done

. /etc/profile.d/z00_lmod.sh
module purge
#module load  intel/19.0.0 intelmpi/19.0.0 hdf5/1.8.21

#module load  intel/19.0.0 intelmpi/19.0.0 hdf5/1.8.21

module load  intel/18.0.1 intelmpi/18.0.1 hdf5/1.8.21 python/3.9.1



python3 "$isolation_root/../workflow_isolation.py" prepare-legacy \
    --root "$isolation_root/.." --directory "$PWD"

rerun=""
mpirun ${PACKHOME}/soft/gizmowjl_MFM_gas_Nstar100/GIZMO  GZWJL.PARAM ${rerun}
#mpirun ${PACKHOME}/soft/gizmowjl_MFM_box_AGS/GIZMO  GZWJL.PARAM ${rerun}
#mpirun ${PACKHOME}/soft/gizmowjl_MFM_AGS_POT/GIZMO  GZWJL.PARAM ${rerun}
#mpirun ${PACKHOME}/soft/gizmowjl_MFM/GIZMO  GZWJL.PARAM ${rerun}
#mpirun /obs/hwang/HFW/gizmowjl_MFM_AGS_NOenergy/GIZMO  GZWJL.PARAM ${rerun}
#mpirun ./gizmowjl_MFM_AnalyticHalo_v3y/GIZMO GZWJL.PARAM ${rerun}
# prolem : mpirun ./gizmowjl_MFM_AnalyticHalo_v3y/GIZMO GZWJL.PARAM ${rerun}

