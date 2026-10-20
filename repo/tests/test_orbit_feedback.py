"""Evolved-feedback regressions; free-flight fixtures are not N-body validation."""
from contextlib import redirect_stdout
from copy import deepcopy
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

SOURCE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SOURCE))
from orbit_calibration_v1 import cli, feedback as fb
from orbit_calibration_v1 import solve_orbit_candidates as solver

KPC_PER_KMS_GYR = 1.022712165045695


def configuration():
    c = json.loads((SOURCE/'orbit_calibration_v1/config.frame_check.example.json').read_text())
    cli.verify_observation_source(SOURCE, c)
    c['feedback'] = dict(comparison_time_Gyr=.2, residual_scales=[.3, .3, 2., 2., .01, .01],
        acceptance_tolerances=[.1, .1, 1., 1., .005, .005], initial_steps=[.5, .5, .5, 1., 1., 1.],
        damping=.0001, trust_radius=1., freeze_physics_inputs=True)
    return c


def history(config, initial=None, run_id='baseline', model='Fornax5026', constant_final=False):
    r, v = np.array(config['target']['r_grf_native_kpc']), np.array(config['target']['v_grf_native_kms'])
    if initial is None:
        initial = np.r_[r-KPC_PER_KMS_GYR*.2*v, v]+[.5, .4, -.3, 2., -1.5, .8]
    initial = np.array(initial, float)
    times = [0., .1, .2]
    frames = []
    for number, t in enumerate(times):
        position = initial[:3]+KPC_PER_KMS_GYR*t*initial[3:]
        velocity = initial[3:]
        if constant_final and number:
            position, velocity = r.copy()+1, v.copy()+1
        frames.append(dict(snapshot=number, time_Gyr=t, r=initial[:3] if number==0 else position,
            v=velocity, state=dict(host=dict(host_center_input_code=[0,0,0], host_bulk_input_code=[0,0,0],
                                            length_to_kpc=1., velocity_to_kms=1.))))
    galaxy = np.r_[initial, np.zeros(6)]
    context = dict(model=model, fixed_input_sha256={k:'a'*64 for k in ['mw_ic','dwarf_ic','cube']})
    return dict(model=model, run_id=run_id, context=context, context_sha256=fb.digest(context), frames=frames,
        initial=initial, inputs={'mergegalaxy.par':[['mw.ini','dwarf.ini'], galaxy.tolist()],
                                  'mergecube.par':[['cube.ini','merged.dat'], np.zeros(12).tolist()]},
        index_path='/synthetic-not-a-real-simulation/index.json', index_sha256=fb.digest(run_id), sources=[],
        physical_to_icrs=config['physical_to_icrs'], verification_scope='SYNTHETIC free-flight fixture')


def response_pairs(c, baseline, constant_final=False):
    result=[]
    for axis in range(6):
        pair=dict(axis=axis)
        for sign, key in [(-1,'minus'),(1,'plus')]:
            p=baseline['initial'].copy(); p[axis]+=sign*c['feedback']['initial_steps'][axis]
            pair[key]=history(c,p,run_id=f'axis{axis}_{key}',constant_final=constant_final)
        result.append(pair)
    return result


def disk_index(directory, c):
    """Small frozen extraction in the real schema, explicitly synthetic."""
    h=history(c)
    inputs=h['inputs']
    case=dict(run_id=h['run_id'], orbit_provenance=dict(parameters_12=inputs['mergegalaxy.par'][1]))
    text={'mergegalaxy.par':'mw.ini dwarf.ini\n'+' '.join(map(str,inputs['mergegalaxy.par'][1]))+'\n',
          'mergecube.par':'cube.ini merged.dat\n'+' '.join(map(str,inputs['mergecube.par'][1]))+'\n',
          'GZWJL.PARAM':'TimeMax 1.0\nInitCondFile initial.IC\n', 'case.json':json.dumps(case)}
    inputs['files']={key:dict(text=value,sha256=__import__('hashlib').sha256(value.encode()).hexdigest()) for key,value in text.items()}
    inputs['fixed_input_sha256']=h['context']['fixed_input_sha256']
    refs={}
    for frame in h['frames']:
        n=frame['snapshot']; cache=directory/f'stars{n}.npz'; sp=directory/f'state{n}.json'
        np.savez(cache,id=np.array([2**63+1,2**63+2],np.uint64),is_old=np.array([True,True]),
                 mass_msun=np.ones(2),dr_native_kpc=np.array([[.1,0,0],[-.1,0,0]]),dv_native_kms=np.zeros((2,3)),
                 r_center_old_kpc=frame['r'],v_center_old_kms=frame['v'])
        state=dict(schema_id='canonical_orbit_state_v1',state='PASS',model=h['model'],numsp=n,
            canonical_observer_id=cli.OBSERVER_ID,observer_config_sha256=cli.OBSERVER_SHA,
            host=dict(canonical_observer_id=cli.OBSERVER_ID,observer_config_sha256=cli.OBSERVER_SHA),
            source_configuration_sha256={key:inputs['files'][key]['sha256'] for key in text if key!='case.json'},
            raw_original_snapshot_sha256=fb.digest(['SYNTHETIC',n]),cache_sha256=cli.sha(cache),
            cohort_file_sha256='b'*64,physical_time_Gyr=frame['time_Gyr'],
            extractor_sha256='e'*64,source_module_sha256={'canonical_observer.py':'f'*64},
            centers=dict(old_full=dict(r_grf_kpc=frame['r'].tolist(),v_grf_kms=frame['v'].tolist(),N=2)))
        sp.write_text(json.dumps(state))
        refs[str(n)]=dict(state_json=sp.name,state_sha256=cli.sha(sp),cache_npz=cache.name)
    index=dict(config={**c,'cases':[dict(model=h['model'],initial=0,latest=2)]},
        server_inputs=dict(models={h['model']:inputs}),states={h['model']:refs})
    path=directory/'index.json'; path.write_text(json.dumps(index))
    return path


class FeedbackTests(unittest.TestCase):
    def setUp(self): self.c=configuration()

    def test_target_point_recovers_all_six_zero_residuals(self):
        target=self.c['target']
        result=fb.evaluate_point(target['r_grf_native_kpc'],target['v_grf_native_kms'],self.c)
        np.testing.assert_allclose(result['residual_vector'],0,atol=1e-9)
        self.assertFalse(result['final_state_registration_applied'])

    def test_sky_residual_wraps_at_zero_ra(self):
        c=json.loads((SOURCE/'orbit_calibration_v1/config.frame_check.example.json').read_text())
        c['target'].update(ra_deg=359.99,dec_deg=0.)
        cli.verify_observation_source(SOURCE,c)
        a=np.radians(.01)
        r=cli.SUN+139.4*np.array(c['physical_to_icrs']).T@np.array([np.cos(a),np.sin(a),0])
        result=fb.evaluate_point(r,c['target']['v_grf_native_kms'],c)
        self.assertAlmostEqual(result['residual_vector'][0],.02,places=9)
        self.assertAlmostEqual(result['sky_separation_deg'],.02,places=9)

    def test_rigid_preview_preserves_measured_speed_and_does_not_set_target_bulk(self):
        _, r, v=solver.rigid_preview_state([50,0,0],[-40,60,0],[0,50,0],[-30,-10,0])
        np.testing.assert_allclose(r,[0,50,0],atol=1e-12)
        self.assertAlmostEqual(np.linalg.norm(v),np.linalg.norm([-40,60,0]))
        self.assertGreater(np.linalg.norm(v-[-30,-10,0]),10)

    def test_one_common_time_and_no_extrapolation(self):
        h=history(self.c)
        self.c['feedback']['comparison_time_Gyr']=.3
        result=fb.evaluate_history(h,self.c)
        self.assertEqual(result['status'],'INSUFFICIENT_TIME_COVERAGE')
        self.assertNotIn('residual_vector',result)
        self.assertFalse(result['evolved_6D_match_validated'])

    def test_nonadjacent_interpolation_is_rejected(self):
        h=history(self.c); h['frames'].pop(1)
        with self.assertRaisesRegex(ValueError,'adjacent'): fb.state_at(h,.1)

    def test_missing_tolerances_never_certifies_a_match(self):
        del self.c['feedback']['acceptance_tolerances']
        result=fb.evaluate_history(history(self.c),self.c)
        self.assertEqual(result['status'],'MEASURED_NO_ACCEPTANCE_TOLERANCES')
        self.assertFalse(result['evolved_6D_match_validated'])

    def test_missing_fixed_input_hashes_refuses_response_plan(self):
        h=history(self.c); h['context']['fixed_input_sha256']={}
        with self.assertRaisesRegex(ValueError,'fixed IC inputs'): fb.plan_feedback(h,self.c)

    def test_plan_does_not_mutate_history_or_launch_simulation(self):
        h=history(self.c); before=h['initial'].copy()
        result=fb.plan_feedback(h,self.c)
        self.assertEqual(len(result['trials']),12)
        self.assertFalse(result['simulations_submitted'])
        np.testing.assert_array_equal(h['initial'],before)
        for p in result['trials']:
            np.testing.assert_array_equal(p['merger']['parameters_12'][6:],np.zeros(6))

    def test_actual_free_flight_responses_improve_forward_residual(self):
        h=history(self.c); result=fb.propose_update(h,response_pairs(self.c,h),self.c)
        actual=history(self.c,result['proposed_initial'],run_id='verified-new-evolution')
        measurement=fb.evaluate_history(actual,self.c)
        verified=fb.verify_proposal(actual,measurement,result,self.c)
        self.assertTrue(verified['accept_for_next_iteration'])
        self.assertLess(verified['actual_score'],verified['previous_score'])
        self.assertFalse(result['evolved_6D_match_validated'])
        self.assertEqual(result['status'],'AWAITING_FORWARD_EVOLUTION')

    def test_five_pairs_or_duplicate_axis_cannot_fake_complete_response(self):
        h=history(self.c); pairs=response_pairs(self.c,h)
        for wrong in [pairs[:5],pairs[:5]+[pairs[0]]]:
            with self.assertRaisesRegex(ValueError,'Six distinct'): fb.propose_update(h,wrong,self.c)

    def test_cross_model_trials_and_reused_run_ids_are_rejected(self):
        h=history(self.c)
        for field,value in [('model','Fornax7004'),('run_id','baseline')]:
            pairs=response_pairs(self.c,h); pairs[0]['plus'][field]=value
            with self.assertRaises(ValueError): fb.propose_update(h,pairs,self.c)

    def test_rank_deficient_evolved_response_is_rejected(self):
        h=history(self.c,constant_final=True)
        with self.assertRaisesRegex(ValueError,'cannot constrain'): fb.propose_update(h,response_pairs(self.c,h,True),self.c)

    def test_verification_requires_new_run_with_proposed_initial_state(self):
        h=history(self.c); result=fb.propose_update(h,response_pairs(self.c,h),self.c)
        with self.assertRaises(ValueError): fb.verify_proposal(h,fb.evaluate_history(h,self.c),result,self.c)
        other=history(self.c,run_id='wrong-new-run')
        with self.assertRaises(ValueError): fb.verify_proposal(other,fb.evaluate_history(other,self.c),result,self.c)

    def test_worse_actual_evolution_is_rejected(self):
        h=history(self.c); p=fb.propose_update(h,response_pairs(self.c,h),self.c)
        actual=history(self.c,p['proposed_initial'],run_id='bad-new-evolution')
        actual['frames'][-1]['v']=actual['frames'][-1]['v']+30
        measured=fb.evaluate_history(actual,self.c)
        result=fb.verify_proposal(actual,measured,p,self.c)
        self.assertFalse(result['accept_for_next_iteration'])
        self.assertEqual(result['status'],'REJECTED_NO_ACTUAL_IMPROVEMENT')

    def test_frozen_loader_cli_and_output_no_overwrite(self):
        with tempfile.TemporaryDirectory() as name:
            directory=Path(name); index=disk_index(directory,self.c)
            h=fb.load_history(index,'Fornax5026')
            self.assertEqual(len(h['frames']),3)
            config=directory/'config.json'; config.write_text(json.dumps(self.c))
            args=['--config',str(config),'--stage','evaluate','--index',str(index),
                  '--model','Fornax5026','--run-directory',str(directory/'out')]
            with redirect_stdout(io.StringIO()): self.assertEqual(cli.main(args),0)
            result=json.loads((directory/'out/FEEDBACK_RESULT.json').read_text())
            self.assertEqual(result['status'],'EVOLVED_6D_MISMATCH')
            with self.assertRaises(FileExistsError),redirect_stdout(io.StringIO()): cli.main(args)
            # Cache content change must be refused rather than silently trusted.
            with (directory/'stars1.npz').open('ab') as stream: stream.write(b'changed')
            with self.assertRaisesRegex(ValueError,'changed native cache'): fb.load_history(index,'Fornax5026')

    def test_extract_can_start_before_a_target_crossing_exists(self):
        c=deepcopy(self.c); c['cases']=[dict(model='Fornax7004',initial=0,latest=59,extract_snapshots=[0,58,59])]
        result=cli.validate_config(c,require_brackets=False)
        self.assertEqual(result['cases'][0]['extract_snapshots'],[0,58,59])
        with self.assertRaises(ValueError): cli.validate_config(c)

    def test_joint_time_search_uses_one_epoch_and_reports_boundary(self):
        h=history(self.c)
        del self.c['feedback']['comparison_time_Gyr']
        self.c['feedback']['search_time_window_Gyr']=[.1,.2]
        result=fb.evaluate_history(h,self.c)
        self.assertTrue(.1<=result['comparison_time_Gyr']<=.2)
        self.assertIn('Joint six-residual',result['time_selection'])
        with self.assertRaisesRegex(ValueError,'Freeze one'): fb.plan_feedback(h,self.c)

    def test_partial_time_window_never_claims_best_epoch(self):
        h=history(self.c)
        del self.c['feedback']['comparison_time_Gyr']
        self.c['feedback']['search_time_window_Gyr']=[.1,.3]
        result=fb.evaluate_history(h,self.c)
        self.assertEqual(result['status'],'INSUFFICIENT_TIME_COVERAGE')
        self.assertNotIn('comparison_time_Gyr',result)

    def test_editing_state_centre_without_matching_frozen_sha_is_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            directory=Path(name); index=disk_index(directory,self.c)
            p=directory/'state1.json'; state=json.loads(p.read_text())
            state['centers']['old_full']['v_grf_kms'][0]+=10
            p.write_text(json.dumps(state))
            with self.assertRaisesRegex(ValueError,'State SHA'): fb.load_history(index,'Fornax5026')

    def test_legacy_identity_is_bound_to_actual_initial_raw_hash(self):
        with tempfile.TemporaryDirectory() as name:
            directory=Path(name); path=disk_index(directory,self.c); index=json.loads(path.read_text())
            item=index['server_inputs']['models']['Fornax5026']['files']['case.json']
            case=json.loads(item['text']); case.update(identity_kind='content_addressed_legacy_run',initial_raw_snapshot_sha256='0'*64)
            item['text']=json.dumps(case); item['sha256']=__import__('hashlib').sha256(item['text'].encode()).hexdigest()
            path.write_text(json.dumps(index))
            with self.assertRaisesRegex(ValueError,'actual snapshot 0'): fb.load_history(path,'Fornax5026')

    def test_capture_fixed_input_hashes_and_explicit_analysis_identity_without_model_edit(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name); model=root/'Fornax5026'; model.mkdir()
            for key in ['mw_ic','dwarf_ic','cube']:
                (root/(key+'.ini')).write_text('fixed input '+key)
            identity=dict(run_id='explicit-legacy',**{k:str(root/(k+'.ini')) for k in ['mw_ic','dwarf_ic','cube']})
            sidecar=root/'analysis-run-identity.json'; sidecar.write_text(json.dumps(identity))
            for name in ['mergegalaxy.par','mergecube.par']:
                (model/name).write_text('first.ini second.ini\n'+' '.join(['0']*12)+'\n')
            (model/'GZWJL.PARAM').write_text('TimeMax 3.3\n')
            record=cli.capture_inputs(root,[dict(model='Fornax5026',run_identity_file=str(sidecar))],True)['models']['Fornax5026']
            self.assertEqual(set(record['fixed_input_sha256']),{'mw_ic','dwarf_ic','cube'})
            self.assertFalse((model/'case.json').exists())
            self.assertEqual(record['files']['case.json']['text'],sidecar.read_text())


if __name__=='__main__': unittest.main()
