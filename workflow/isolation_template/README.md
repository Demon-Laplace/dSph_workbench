# Isolation shell templates

Deploy these five launchers into the server's unnumbered `isoFornax/Fornax/`
template. Keep its existing `GZWJL.PARAM`. Python implementations live in
`isoFornax/python/` and the data-root `workflow_isolation.py`; do not copy them
into new model directories.

For manual preparation, create `isoFornax/FornaxNUMBER/`, copy GZWJL.PARAM and
only the needed launchers, and upload `FornaxNUMBER.ini`. Submit subastrothin.sh
or submeso1.sh from that model directory. These preserve the historical Slurm
resources/solver and first call prepare-legacy to validate the numbered input,
materialize TREECOOL and convert the IC. Existing output is never restarted by
this preparation command. Configure the converter and cooling template using
the existing isolation_profiles.default server setting.

Prefer the local workflow.py isolation MODEL [MODEL ...] --name RUN --submit
command for model/parameter/hash verification and automatic upload. It creates
separate frozen runs rather than modifying these manual model directories.
