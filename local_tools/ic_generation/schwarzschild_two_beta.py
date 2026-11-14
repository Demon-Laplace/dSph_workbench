"""Compatibility reader for one stellar component with two beta windows.

Imports the original solver without editing it. This file does nothing until
explicitly run. It never submits jobs, publishes, or starts an evolution run.
"""
import argparse,ast,configparser,hashlib,importlib.util,json,os,re,tempfile
from pathlib import Path
import numpy as np

REGIONS=('outer','inner')

def grid(value,agama_module):
    value=value.strip()
    if '(' in value:
        node=ast.parse(value,mode='eval').body
        if not isinstance(node,ast.Call) or node.keywords:
            raise ValueError('Use a literal grid or an allowed numeric grid function.')
        name=ast.unparse(node.func)
        allowed={'numpy.linspace':np.linspace,'np.linspace':np.linspace,
            'numpy.geomspace':np.geomspace,'np.geomspace':np.geomspace,
            'agama.nonuniformGrid':agama_module.nonuniformGrid}
        if name not in allowed:raise ValueError('Unsupported grid function: '+name)
        values=allowed[name](*[ast.literal_eval(a) for a in node.args])
    else:
        values=[float(x) for x in re.split(r'[,\s]+',value) if x]
    values=np.asarray(values,dtype=float)
    if values.ndim!=1 or len(values)<2 or not np.isfinite(values).all() or np.any(np.diff(values)<=0) or values[0]<0:
        raise ValueError('Grid must be finite, nonnegative and strictly increasing.')
    return values

def prepare_config(path,agama_module):
    config=configparser.RawConfigParser()
    if not config.read(path):raise FileNotFoundError(path)
    specs={}
    for section in config.sections():
        if section.lower().startswith('component'):
            value=config[section]
            dual=value.get('beta','').strip().lower()=='two_zone'
            has_fields=any(k in value for k in ('beta_inner','beta_outer','kinemgrid_inner','kinemgrid_outer'))
            if has_fields and not dual:raise ValueError('Dual-beta fields require beta=two_zone.')
            if dual:
                if 'kinemgrid' in value:raise ValueError('Remove the old single kinemgrid in two_zone mode.')
                entry=[]
                for region in REGIONS:
                    beta=float(value['beta_'+region]);radius=grid(value['kinemgrid_'+region],agama_module)
                    if not np.isfinite(beta) or beta>=1:raise ValueError('beta must be finite and below 1.')
                    entry.append(dict(region=region,beta=beta,grid=radius,
                        degree=int(value.get('kinemdegree','2')),normalize=(region=='inner')))
                if entry[1]['grid'][-1]>=entry[0]['grid'][0]:
                    raise ValueError('Inner and outer supports must leave a non-overlapping transition.')
                specs[section.lower()]=entry
                config.remove_option(section,'beta')
                for region in REGIONS:
                    config.remove_option(section,'beta_'+region)
                    config.remove_option(section,'kinemgrid_'+region)
                config.set(section,'kinemgrid',', '.join(f'{v:.16g}' for v in entry[0]['grid']))
            for key in ('gridr','gridz','kinemgrid'):
                if key in value:
                    config.set(section,key,', '.join(f'{v:.16g}' for v in grid(value[key],agama_module)))
    if not specs:raise ValueError('No beta=two_zone component was specified.')
    return config,specs

def load_legacy(path):
    spec=importlib.util.spec_from_file_location('_original_schwarzschild',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def createModel(ini_file,base_solver=None):
    import agama
    if base_solver is None:
        candidate=Path(__file__).resolve().with_name('schwarzschild.py')
        base_solver=candidate if candidate.is_file() else None
    if base_solver is None:raise FileNotFoundError('Provide --base-solver or place this file beside schwarzschild.py.')
    legacy=load_legacy(base_solver);config,specs=prepare_config(ini_file,agama)
    with tempfile.NamedTemporaryFile(mode='w',suffix='.ini',prefix='two_beta_',dir=Path.cwd(),delete=False) as f:
        config.write(f);temporary=Path(f.name)
    try:model=legacy.createModel(str(temporary))
    finally:temporary.unlink()
    for name,regions in specs.items():
        comp=model.components[name]
        comp.targets=[comp.targets[0]]+[agama.Target(type='KinemShell',gridr=r['grid'],degree=r['degree']) for r in regions]
        comp.beta_regions=regions
        if hasattr(comp,'beta'):del comp.beta
    model.original_solver=legacy
    return model

def moment_constraint(matrix,mass,beta,normalize):
    matrix=np.asarray(matrix,dtype=np.float64)
    if (matrix.ndim!=2 or matrix.shape[1]==0 or matrix.shape[1]%2
            or not np.isfinite(matrix).all() or not np.isfinite(mass) or mass<=0):
        raise ValueError('Finite radial/tangential moment pairs and positive mass required.')
    n=matrix.shape[1]//2
    radial=matrix[:,:n];tangent=matrix[:,n:]
    equation=2*(1-beta)*radial.T-tangent.T
    if normalize:
        scale=mass*np.mean(2*np.abs(radial)+np.abs(tangent),axis=0)
        if np.any(scale<=0):raise ValueError('Orbit library has an empty kinetic basis.')
        equation=equation/scale[:,None]
    return equation

def runComponent(comp,potential,legacy):
    import agama
    if not hasattr(comp,'beta_regions'):
        return legacy.runComponent(comp,potential)
    kwargs=dict(potential=potential,ic=comp.ic,time=comp.inttime,Omega=comp.Omega,targets=comp.targets)
    result=agama.orbit(**kwargs,**({'dtype':object} if hasattr(comp,'nbody') else {}))
    mass=comp.density.totalMass();density=result[0];tau=comp.targets[0](comp.density)
    matrices=[density.T];rhs=[tau];penalty=[np.ones(len(tau))/(mass/len(tau))]
    for index,region in enumerate(comp.beta_regions,1):
        equation=moment_constraint(result[index],mass,region['beta'],region['normalize'])
        matrices.append(equation);rhs.append(np.zeros(equation.shape[0]));penalty.append(np.ones(equation.shape[0]))
    matrices.append(np.ones((1,len(comp.ic))));rhs.append(np.array([mass]));penalty.append(np.array([np.inf]))
    regularizer=np.ones(len(comp.ic))/(mass/len(comp.ic))**2/len(comp.ic)*.1
    weights=agama.solveOpt(matrix=matrices,rhs=rhs,rpenl=penalty,xpenq=regularizer)
    weights=np.asarray(weights,dtype=float)
    if weights.shape!=(len(comp.ic),) or not np.isfinite(weights).all():
        raise RuntimeError('Orbit solver returned invalid weights.')
    diagnostics=dict(density_half_L1=float(.5*np.sum(abs(weights@density-tau))/mass),
        relative_mass_error=float(abs(weights.sum()/mass-1)),minimum_weight=float(weights.min()),regions=[])
    for index,region in enumerate(comp.beta_regions,1):
        physical=np.asarray(weights@result[index],float);n=len(physical)//2
        if np.any(physical[:n]<=0):raise RuntimeError('No positive radial moment in a beta basis.')
        actual=1-physical[n:]/(2*physical[:n]);error=float(np.max(abs(actual-region['beta'])))
        diagnostics['regions'].append(dict(region=region['region'],target_beta=region['beta'],
            grid_kpc=region['grid'].tolist(),actual_beta_per_basis=actual.tolist(),maximum_error=error))
    good=(diagnostics['density_half_L1']<.002 and diagnostics['relative_mass_error']<1e-6
        and diagnostics['minimum_weight']>=0 and all(r['maximum_error']<.03 for r in diagnostics['regions']))
    comp.two_beta_validation=diagnostics
    if not good:raise RuntimeError('Two-beta constraints did not pass: '+json.dumps(diagnostics))
    if hasattr(comp,'nbody'):
        status,particles=agama.sampleOrbitLibrary(comp.nbody,result[-1],weights)
        if not status:raise RuntimeError('Failed to sample the orbit library.')
        comp.nbodymodel=particles
    comp.weights=weights;comp.densitydata=density;comp.kinemdata=result[1];comp.innerkinemdata=result[2]

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('ini',type=Path)
    parser.add_argument('--base-solver',type=Path);parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--seed',type=int,default=42);args=parser.parse_args()
    source=args.ini.resolve();base=args.base_solver.resolve() if args.base_solver else None
    if not source.is_file():parser.error('INI source must exist.')
    if base is not None and not base.is_file():parser.error('--base-solver must exist.')
    if args.output.exists():parser.error('--output must be a new directory.')
    if not 0<=args.seed<2**31:parser.error('--seed must be in [0, 2**31).')
    import agama
    args.output.mkdir(parents=True)
    os.chdir(args.output);np.random.seed(args.seed);agama.setRandomSeed(args.seed)
    source_bytes=source.read_bytes();Path('parameters.ini').write_bytes(source_bytes)
    receipt=dict(state='constructing',seed=args.seed,source_ini=str(source),
        source_ini_sha256=hashlib.sha256(source_bytes).hexdigest(),
        compatibility_reader_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        no_evolution_simulation=True)
    Path('two_beta_build.json').write_text(json.dumps(receipt,indent=2))
    try:
        model=createModel(source,base)
        original_path=Path(model.original_solver.__file__).resolve()
        receipt.update(original_solver=str(original_path),
            original_solver_sha256=hashlib.sha256(original_path.read_bytes()).hexdigest())
        Path('two_beta_build.json').write_text(json.dumps(receipt,indent=2))
        for name,comp in model.components.items():
            try:
                runComponent(comp,model.potential,model.original_solver)
            finally:
                if hasattr(comp,'two_beta_validation'):
                    Path('model_'+name+'_two_beta_validation.json').write_text(json.dumps(comp.two_beta_validation,indent=2))
            if hasattr(comp,'nbody'):agama.writeSnapshot('model_'+name+'.nbody',comp.nbodymodel)
            values={k:getattr(comp,k) for k in ('ic','inttime','weights','densitydata','kinemdata','innerkinemdata') if hasattr(comp,k)}
            np.savez_compressed('model_'+name+'.data',**values)
    except Exception as error:
        receipt.update(state='failed',error=str(error))
        Path('two_beta_build.json').write_text(json.dumps(receipt,indent=2))
        raise
    receipt['state']='complete';Path('two_beta_build.json').write_text(json.dumps(receipt,indent=2))
    print('Static particle construction complete. No evolution has been started.')

if __name__=='__main__':main()
