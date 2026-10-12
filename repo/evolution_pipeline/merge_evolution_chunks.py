#!/usr/bin/env python3
"""Merge independently extracted snapshot chunks and reuse the main derivation path."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pipeline-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--chunk", type=Path, action="append", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sys.path.insert(0, str(args.pipeline_root.resolve()))
    import fornax_evolution_pipeline as pipeline

    config, config_path = pipeline.load_config(args.config)
    paths = pipeline.resolve_paths(config, config_path)
    sources = [paths["csv"], *(path.resolve() for path in args.chunk)]
    for source in sources:
        if not source.exists():
            raise FileNotFoundError(source)

    frames = [pd.read_csv(source) for source in sources]
    combined = pd.concat(frames, ignore_index=True, sort=False)
    membership_counts = set(
        combined["cgm_dwarf_gas_membership_count"].dropna().astype(int)
    )
    membership_hashes = set(
        combined["cgm_dwarf_gas_membership_sha256"].dropna().astype(str)
    )
    if len(membership_counts) != 1 or len(membership_hashes) != 1:
        raise RuntimeError(
            "Chunk CGM-membership definitions differ: "
            f"counts={membership_counts}, hashes={membership_hashes}"
        )

    overlap_columns = [
        "time_gyr",
        "dwarf_center_x_kpc",
        "dwarf_center_y_kpc",
        "dwarf_center_z_kpc",
        "gas_mass_msun",
        "hi_mass_msun",
        "re_major_kpc",
        "sigma_los_kms",
        "cgm_density_g_cm3",
        "v_rel_cgm_kms",
        "ram_pressure_dyn_cm2",
    ]
    for snapshot, duplicate in combined.loc[
        combined.duplicated("snapshot", keep=False)
    ].groupby("snapshot"):
        values = duplicate[overlap_columns].to_numpy(dtype=float)
        if not np.allclose(values, values[0], rtol=1.0e-11, atol=1.0e-13, equal_nan=True):
            raise RuntimeError(f"Non-identical duplicate measurements for snapshot {snapshot}")

    merged = (
        combined.sort_values(["snapshot"])
        .drop_duplicates("snapshot", keep="last")
        .reset_index(drop=True)
    )
    expected = [number for number, _ in pipeline.discover_snapshots(paths["snapshot_dir"], config)]
    actual = merged["snapshot"].astype(int).tolist()
    if actual != expected:
        missing = sorted(set(expected) - set(actual))
        extra = sorted(set(actual) - set(expected))
        raise RuntimeError(f"Incomplete merged sequence: missing={missing}, extra={extra}")

    merged["analysis_config_sha256"] = pipeline.analysis_config_hash(config)
    merged = pipeline.add_derived_columns(merged, config)
    pipeline.atomic_write_csv(merged, paths["csv"])
    pipeline.save_metadata(merged, config, config_path, paths)
    outputs = pipeline.plot_timeseries(config, config_path)
    print(
        f"merged_rows={len(merged)} membership_count={next(iter(membership_counts))}",
        flush=True,
    )
    for output in outputs:
        print(output, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
