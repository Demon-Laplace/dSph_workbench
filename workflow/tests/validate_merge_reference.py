#!/usr/bin/env python3
"""Compare one-pass IC with the actual legacy two-pass pipeline in scratch space.
Run on a compute node for production-size inputs. Never uses a model as scratch.
Deletes only its own generated ASCII/binaries after complete equality validation.
"""
import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys
import time
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import workflow_merge as merge
import workflow_runtime as runtime


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--galaxy',type=Path,required=True); p.add_argument('--cube',type=Path,required=True)
    p.add_argument('--merger',type=Path,required=True); p.add_argument('--converter',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True)
    a=p.parse_args(); work=a.work.resolve()
    if work.exists(): raise ValueError('Choose a new scratch directory')
    g,c=merge.read_par(a.galaxy),merge.read_par(a.cube)
    if Path(c[0][1]) not in (Path('mergegalaxy.dat'),Path('mergegalaxye.dat')): raise ValueError('Unexpected merge chain')
    inputs=[(a.galaxy.parent/g[0][0]).resolve(),(a.galaxy.parent/g[0][1]).resolve(),(a.cube.parent/c[0][0]).resolve()]
    if not all(x.is_file() for x in inputs): raise ValueError('Missing inputs')
    work.mkdir(parents=True)
    for name,source in zip(('mw.ini','dwarf.ini','cube.ini'),inputs): (work/name).symlink_to(source)
    (work/'galaxy.par').write_text('mw.ini dwarf.ini\n'+' '.join(g[1])+'\n')
    (work/'cube.par').write_text('cube.ini galaxy.dat\n'+' '.join(c[1])+'\n')
    source=Path(merge.__file__).with_name('ic_merge.c'); executable=work/'ic_merge'
    subprocess.run(['gcc','-O2','-std=c99','-ffp-contract=off','-fno-fast-math',str(source),'-lm','-o',str(executable)],check=True)
    report={'started_utc':runtime.now(),'inputs':{str(f):merge.sha(f) for f in inputs},
            'galaxy_parameters':g,'cube_parameters':c,
            'merger_sha256':merge.sha(a.merger),'converter_sha256':merge.sha(a.converter),'source_sha256':merge.sha(source),'seconds':{}}
    def run(label,cmd):
        start=time.monotonic()
        with (work/(label+'.log')).open('w') as log: subprocess.run(list(map(str,cmd)),cwd=work,stdout=log,stderr=subprocess.STDOUT,check=True)
        report['seconds'][label]=time.monotonic()-start
    run('old_galaxy',[a.merger,'galaxy.par','galaxy.dat'])
    run('old_cube',[a.merger,'cube.par','old.dat'])
    run('onepass',[executable,'galaxy.par','cube.par','new.dat'])
    report['ascii_sha256']={name:merge.sha(work/name) for name in ('old.dat','new.dat')}
    if len(set(report['ascii_sha256'].values()))!=1: raise ValueError('ASCII differs from legacy result')
    for prefix in ('old','new'): run(prefix+'_convert',[a.converter,prefix+'.dat',prefix+'.IC_000'])
    report['binary_sha256']={name:merge.sha(work/name) for name in ('old.IC_000','new.IC_000')}
    if len(set(report['binary_sha256'].values()))!=1: raise ValueError('Binary differs from legacy result')
    report['validation']=merge.validate_binary(work/'new.dat',work/'new.IC_000')
    report['byte_identical_ascii_and_binary']=True
    report['finished_utc']=runtime.now()
    runtime.write_json(work/'validation.json',report)
    for name in ('galaxy.dat','old.dat','new.dat','old.IC_000','new.IC_000'): (work/name).unlink()
    report['temporary_particle_outputs_deleted']=True
    runtime.write_json(work/'validation.json',report)
    print(json.dumps(report,indent=2))


if __name__=='__main__': main()
