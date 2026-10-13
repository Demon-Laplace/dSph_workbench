# Shared legacy isolation Python

These 14 files are exact copies of the identical Python files found in 84
valid numerical server isolation models on 2026-10-13. No numerical changes
were made. The server deployment lives at `data/isoFornax/python/`; each model
keeps relative symlinks to preserve its existing shell commands.

Keep this deployment depth. `check_density_profile.py` derives `isoFornax/`
and `data/` using its resolved source path and imports existing helpers from
`isoFornax/` first. Those root helpers can differ from this historical set and
must not be silently replaced. Other legacy entry points use the shared set.
Run from the model working directory so inputs and products remain model-local.

Edits to this shared code affect every linked legacy model and require explicit
scientific validation. Frozen new workflow runs and the historical template
are not migrated by this change. Never edit a user's existing `.sh` to adapt
it to this layout without an explicit request.
