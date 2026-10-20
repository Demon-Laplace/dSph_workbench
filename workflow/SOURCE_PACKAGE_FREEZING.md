# Canonical orbit package source handling

`orbit_calibration_v1/` is an explicit production-source package. Source export
and analysis `prepare` share its file policy: Python modules, `README` or
`README.md`, and `config*.json`, recursively through ordinary source directories.
Python and JSON must parse, and each allowed file must be at most 2 MiB.

Hidden entries and the `__pycache__`, `cache`, `data`, `figures`, `logs`, `output`,
`outputs`, and `results` directories are excluded. Other data files are excluded.
Links in the selected package are rejected rather than followed. No particle
input, generated result or cached analysis is copied into the code snapshot.

`prepare` copies accepted package files under `run/code/orbit_calibration_v1/`
and records each file in `manifest.json.code_sha256`. The existing input verifier
rejects changes to those frozen files. Existing info/plot output locations and
old frozen tasks retain their behavior; an absent optional package does not
change how older production code is prepared. Git releases already include
tracked nested files, while source export and prepare require this explicit
allowlist. Install updated workflow modules separately from activating a release.

Run the focused synthetic checks from `workflow/`:

```bash
python3 -m unittest discover -s tests -p test_orbit_package_freezing.py
```

On Windows the prepare fixture substitutes copies for its input links when
testing package copying. Source import, manifest checks and tamper detection
remain real; link-rejection cases use actual links and skip only when the host
cannot create them. Linux runs use real prepare input links. Set
`DSPH_TEST_SCRATCH` to an existing writable scratch directory when the host's
default temporary directory is unavailable.
