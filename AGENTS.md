# dSph workbench agent entry point

Read [workflow/AGENT_HANDOFF.md](workflow/AGENT_HANDOFF.md) before operating the
local/IC/server workflow. It records the directory conventions, server write
boundary, startup-retry cleanup policy, density analysis and publication steps.
Use workflow/WORKFLOW.md for further command details. Inspect installed copies
and Git status before changing code. Do not edit frozen simulation runs.

Submit checks for multiple simulations as one Slurm job on one node, even when
their isolation run names differ. Use repeated `--run MODEL=RUN` arguments to
identify those runs. Process models sequentially with the configured worker count
(normally 16), including video-only checks. Do not submit one job per model or
run multiple worker groups concurrently. Preserve already submitted tasks.

Production orbital simulations use astro_fat, one separate job per model,
through the existing subastrofat.sh with its current resource settings. For
Fornax7004/7005, use snap2csv exports of their completed isolation snapshot 50
and retain timescale_fd=16.6. Validate exports and merged ICs before submission;
preserve original snapshots and existing user shell scripts.
Fornax7004/7005 orbital and orientation transforms must match Fornax5028.
Use automatic startup monitoring: after 140 running seconds and nonzero
integration progress, exit the monitor and leave the simulation running.
Avoid repeated manual checks or confirmation requests.
After successful startup detection, return simulation details (IC/model and
Milky Way/background inputs, particle counts/mass resolution, orbit/orientation,
physics parameters and job resources) and end monitoring.
Integrate monitoring into the existing batch allocation as a sleeping child,
without an extra Slurm job or reserved CPU. Use a generated wrapper, preserving
user launchers. Queue time uses no monitor node resources; local disconnection
must not stop server detection or simulation.
