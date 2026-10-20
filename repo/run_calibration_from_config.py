#!/usr/bin/env python3
"""Run the original orbit entry with an explicit canonical six-dimensional configuration.

Legacy response-table configurations are not silently reinterpreted.
"""
from pathlib import Path
import argparse,json,subprocess,sys

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config',type=Path)
    args,forward=parser.parse_known_args()
    path=args.config.resolve()
    config=json.loads(path.read_text(encoding='utf-8-sig'))
    if config.get('schema_id')!='canonical_orbit_recalibration_v1':
        raise SystemExit('Old CSV response coefficients lack a certified observer and H/B origin. Supply a new canonical ICRS configuration; preserved historical configurations are not converted.')
    if config.get('target',{}).get('input_frame')!='icrs_heliocentric':
        raise SystemExit('Declare target.input_frame=icrs_heliocentric. Li Table 2 GC radial velocities are not heliocentric LOS inputs.')
    return subprocess.run([sys.executable,'-B',str(Path(__file__).with_name('calibrate_orbit.py')),'--config',str(path),*forward],check=False).returncode

if __name__=='__main__':raise SystemExit(main())
