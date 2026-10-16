"""Run a frozen isolation case with the configured native MPI solver on this host."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import workflow_runtime as runtime
import workflow_isolation as isolation


def directory(root, model, name):
    runtime.label(model); runtime.label(name)
    return runtime.confined(root, root / 'workspace/local-isolation' / model / name)


def progress(run, status, samples):
    log = run / 'solver.log'
    if not log.exists(): return
    with log.open('rb') as f:
        f.seek(max(0, log.stat().st_size - 1024*1024))
        tail = f.read().decode(errors='replace')
    matches = re.findall(r'Sync-Point \d+, Time:\s*([0-9.eE+-]+)', tail)
    if not matches: return
    current = float(matches[-1]); elapsed = time.time() - status['started_epoch']
    if not samples or current > samples[-1]['simulation_time']:
        samples.append({'elapsed_seconds': elapsed, 'simulation_time': current})
    del samples[:-40]
    value = {'updated_utc': runtime.now(), 'elapsed_seconds': elapsed,
             'simulation_time': current, 'target_time': status['target_time'], 'samples': samples}
    # Exclude initialization and use the recent measured integration rate.
    useful = [p for p in samples if p['simulation_time'] > 0]
    if len(useful) >= 2 and useful[-1]['simulation_time'] > useful[0]['simulation_time']:
        dt = useful[-1]['elapsed_seconds'] - useful[0]['elapsed_seconds']
        speed = (useful[-1]['simulation_time'] - useful[0]['simulation_time']) / dt
        value['estimated_remaining_hours'] = max(0, (status['target_time']-current)/speed/3600)
        value['estimate_scope'] = 'recent integration speed; excludes future speed changes and suspension'
    runtime.write_json(run / 'progress.json', value)


def execute(root, run):
    with (run / 'worker.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        request = json.loads((run/'request.json').read_text())
        status = json.loads((run/'status.json').read_text())
        if status['state'] != 'starting': raise ValueError('Local run is not waiting for execution')
        status.update(state='running', started_utc=runtime.now(), started_epoch=time.time(),
                      pid=os.getpid(), target_time=request['time_max'])
        runtime.write_json(run/'status.json', status)
        environment = dict(os.environ, OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                           NUMEXPR_NUM_THREADS='1', TMPDIR=str(run/'.tmp'))
        (run/'.tmp').mkdir(exist_ok=True)
        try:
            for name, digest in request['frozen_sha256'].items():
                if isolation.sha(run/name) != digest: raise ValueError('Changed frozen file: '+name)
            os.nice(10)
            with (run/'converter.log').open('w') as log:
                subprocess.run([str(run/'IC_ascii2snap'),str(run/'Fornax.ini'),str(run/'Fornax.IC_000')],
                               cwd=run,env=environment,stdout=log,stderr=subprocess.STDOUT,check=True)
            command = [x.format(solver=str(run/'GIZMO'), parameters='GZWJL.PARAM',
                                processes=request['processes']) for x in request['solver_command']]
            with (run/'solver.log').open('w') as log, subprocess.Popen(command,cwd=run,env=environment,
                    stdout=log,stderr=subprocess.STDOUT) as process:
                samples=[]
                while True:
                    try: code=process.wait(timeout=15)
                    except subprocess.TimeoutExpired:
                        progress(run,status,samples); continue
                    progress(run,status,samples)
                    if code: raise subprocess.CalledProcessError(code,command)
                    break
            status.update(state='complete', verification=isolation.verify_completion(run))
        except Exception as error:
            status.update(state='failed',error=str(error),finished_utc=runtime.now())
            runtime.write_json(run/'status.json',status)
            raise
        status['finished_utc']=runtime.now()
        runtime.write_json(run/'status.json',status)


def start(args):
    import workflow_models as models
    models.model_name(args.model)
    cfg=runtime.config(args.root)['isolation_local']
    source=models.find_initial(args,args.ic_run,args.model)
    build=json.loads((source/'build.json').read_text())
    if (build['state']!='complete' or Path(build['source_ini']).stem!='IC_'+args.model or
        isolation.sha(source/'Fornax.ini')!=build['particle_sha256'] or
        isolation.sha(source/'parameters.ini')!=build['ini_sha256']):
        raise ValueError('IC build identity/checksum mismatch')
    processes=args.processes or cfg.get('processes',4)
    if processes < 1 or processes > (os.cpu_count() or 1): raise ValueError('Invalid MPI process count')
    if args.time_max is not None and args.time_max <= 0: raise ValueError('--time-max must be positive')
    run=directory(args.root,args.model,args.name)
    identity={'model':args.model,'ic_run':args.ic_run,'processes':processes,
              'particle_sha256':build['particle_sha256'],'time_max':args.time_max,
              'profile':cfg}
    if args.dry_run:
        print(json.dumps({'directory':str(run),'identity':identity},indent=2));return 0
    run.parent.mkdir(parents=True,exist_ok=True)
    with (run.parent/(run.name+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if run.exists():
            if json.loads((run/'identity.json').read_text())!=identity:
                raise ValueError('Local run name already belongs to different inputs/settings')
            print((run/'status.json').read_text());return 0
        template=(args.root/cfg['template']).resolve()
        files={'Fornax.ini':source/'Fornax.ini','parameters.ini':source/'parameters.ini',
               'GIZMO':(args.root/cfg['solver']).resolve(),'IC_ascii2snap':(args.root/cfg['converter']).resolve(),
               'TREECOOL':template/'TREECOOL'}
        if not all(p.is_file() for p in files.values()):raise ValueError('Missing native solver, converter or template input')
        run.mkdir()
        for name,path in files.items():shutil.copy2(path,run/name)
        for name in ('workflow_local_isolation.py','workflow_runtime.py','workflow_isolation.py'):
            shutil.copy2(Path(__file__).with_name(name),run/name)
        values=dict(cfg.get('parameter_overrides',{}))
        values.update(InitCondFile='Fornax.IC_000',OutputDir='output',ResubmitOn='0',
                      NumFilesWrittenInParallel=str(min(4,2**(processes.bit_length()-1))))
        if args.time_max is not None:
            values['TimeMax']=str(args.time_max)
        isolation.rewrite_parameters(template/'GZWJL.PARAM',run/'GZWJL.PARAM',values)
        fields={line.split('%',1)[0].split()[0]:line.split('%',1)[0].split()[1]
                for line in (run/'GZWJL.PARAM').read_text().splitlines() if len(line.split('%',1)[0].split())>=2}
        (run/'output').mkdir()
        request={'model':args.model,'processes':processes,'time_max':float(fields['TimeMax']),
                 'solver_command':cfg['command'],'frozen_sha256':{p.name:isolation.sha(p) for p in run.iterdir() if p.is_file()}}
        runtime.write_json(run/'request.json',request);runtime.write_json(run/'identity.json',identity)
        runtime.write_json(run/'status.json',{'state':'starting','created_utc':runtime.now(),'directory':str(run)})
        command=[str((args.root/cfg['python']).resolve()),str(run/'workflow_local_isolation.py'),
                 '--root',str(args.root),'--directory',str(run)]
        with (run/'worker.log').open('w') as log:
            process=subprocess.Popen(command,cwd=run,start_new_session=True,stdin=subprocess.DEVNULL,
                                     stdout=log,stderr=subprocess.STDOUT)
        runtime.write_json(run/'process.json',{'pid':process.pid,'started_utc':runtime.now()})
        print(json.dumps({'state':'started','pid':process.pid,'directory':str(run),'processes':processes},indent=2))
    return 0


def status(args):
    run=directory(args.root,args.model,args.name)
    value=json.loads((run/'status.json').read_text())
    if (run/'progress.json').exists():value['progress']=json.loads((run/'progress.json').read_text())
    print(json.dumps(value,indent=2));return 0


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--directory',type=Path,required=True);a=p.parse_args()
    execute(a.root.resolve(),runtime.confined(a.root,a.directory).resolve())
