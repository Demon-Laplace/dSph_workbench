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
