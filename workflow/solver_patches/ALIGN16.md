# Legacy periodic GIZMO: 16-byte allocator

The periodic executable used for Fornax8002 rounds four allocator sizes to eight
bytes, while `drift_particle` uses a 16-byte-aligned SSE store. Changing Slurm's
memory limit cannot repair that mismatch.

`build_align16.py` copies the original source into a new data-workspace directory,
changes only `system/mymalloc.c`, and builds a separately named `GIZMO_align16`:

- Pool base: `posix_memalign(...,16,...)`.
- All four allocation/reallocation paths: round requested sizes to multiples of
  16, with a minimum block size of 16. Subsequent movable-block relocation
  therefore preserves the alignment too.
- Keep the original Config.sh, generated physics macros, compiler and Makefile
  optimization settings. No changes to particle structures, floating-point
  precision, physical equations, ICs, or runtime science parameters.

Run on a compute node with the original Intel MPI/HDF5 module environment:

```bash
python3 build_align16.py --source ORIGINAL_SOURCE \
  --destination NEW_RELEASE/source --jobs 4
```

Do not build inside `/obs/`, replace its original executable, or switch an active
job's executable. The new build needs allocator tests (including relocation),
synthetic startup and binary restart tests. An old-to-new synthetic restart
checks basic checkpoint compatibility; it is not a guarantee for every
production-scale checkpoint or scientific evolution.

## Automatic selection on future production submissions

After the tests succeed, the server's machine-specific `workflow.json` may set:

```json
{
  "production_solver": {
    "executable": "/travail/xhuang/workspace/software/RELEASE/source/GIZMO_align16",
    "build_manifest": "/travail/xhuang/workspace/software/RELEASE/source/aligned-build.json",
    "validation": "/travail/xhuang/workspace/software/RELEASE/validation/validation.json",
    "replaces": [
      "${PACKHOME}/soft/gizmowjl_MFM_box_AGSgas_Nstar100/GIZMO",
      "/obs/xhuang/zimpack/soft/gizmowjl_MFM_box_AGSgas_Nstar100/GIZMO"
    ]
  }
}
```

`simulate --submit` changes the solver path only in the newly frozen launcher's
active command. Existing model/user shells, already submitted jobs, Slurm
resources, `rerun`, parameters and ICs remain intact. The small build/validation
records must declare successful alignment tests before selection. No particle
hashes or routine data certification are performed.

This mapping is deliberately limited to the tested periodic solver family;
isolation executables and unrelated custom solvers are not silently replaced.
Calling `sbatch` directly on an old user shell bypasses workflow selection.

## Validated release (2026-11-21)

Server release: `/travail/xhuang/workspace/software/gizmo-periodic-align16-20261121-v1/`.
Build/test job236891 completed with exit0. Four MPI ranks passed actual allocator
base/allocation/reallocation/relocation tests. Synthetic periodic cases with
4100 particles (unaligned under the old layout) and 4096 particles completed
new-to-new and old-to-new binary restarts from code time1e-5 to2e-5.
Physics macros were unchanged. These are short, single-node tests, not a
96-rank production validation or long-term scientific result.

Before default activation, the current queue was checked: Fornax8002 isolation
job314033 was active; production job314009 remained failed with only its t=0
snapshot. Neither that isolation job nor any submitted production job is
switched or restarted by this deployment.
