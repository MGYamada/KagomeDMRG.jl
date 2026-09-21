#!/usr/bin/env python3
"""Posthoc axial domain diagnostic on the completed N54 comparison.

This is a comparison of local profiles on common interior support, not a
symmetry operation on the open cylinder or a phase/energy selection rule.
"""
import argparse,hashlib,itertools,json,math
from pathlib import Path
import tomllib
root=Path(__file__).resolve().parent.parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('validation',type=Path)
parser.add_argument('output',type=Path)
parser.add_argument('--stage',choices=('release','relax'),default='relax')
args=parser.parse_args()
p=args.validation.resolve()
assert not args.output.exists(), 'refusing existing analysis'

raw=p.read_bytes();d=tomllib.loads(raw.decode());g=d['configuration']
rows={r['branch']:r for r in d['batches'] if r['stage']==args.stage}
assert d['status']=='completed_preparation_comparison_accuracy_separate'
assert d['sources_unchanged'] is True and d['backend_unchanged'] is True
execution=tomllib.loads((p.parent/'execution.toml').read_text())
assert execution['status']=='exited' and execution['worker_exit_confirmed'] and execution['worker_exit_code']==0
assert set(rows)=={'random','hourglass','windmill'}
for r in rows.values():
 assert r['integrity_passed'] is True and r['active_phase']=='completed' and r['lambda']==0 and r['final_nn_model_verified'] is True
 assert r['model_configuration']==g and r['cumulative_sweeps']==(6 if args.stage=='release' else 8)
 for name,sha in r['checkpoint_sha256'].items():assert hashlib.sha256((root/r['checkpoint']/name).read_bytes()).hexdigest()==sha
sites=g['sites'];bonds=g['bonds'];ly=g['Ly'];lookup={(s['x'],s['y'],s['sublattice']):i for i,s in enumerate(sites)}
def bkey(i,j,wy):return (i,j,wy) if i<j else (j,i,-wy)
bl={bkey(b['i']-1,b['j']-1,b['wy']):k for k,b in enumerate(bonds)}
def support(dx):
 ix=[i for i,s in enumerate(sites) if 1<=s['x']<=4 and 1<=s['x']+dx<=4]
 kept=set(ix);ib=[k for k,b in enumerate(bonds) if b['i']-1 in kept and b['j']-1 in kept]
 assert len(ix)==(36 if dx==0 else 27)
 assert len(ib)==(66 if dx==0 else 48)
 return ix,ib

def mapping(dx,dy,ix,ib):
 sm={i:lookup[(sites[i]['x']+dx,(sites[i]['y']+dy)%ly,sites[i]['sublattice'])] for i in ix}
 bm={}
 for k in ib:
  b=bonds[k];i,j=b['i']-1,b['j']-1
  wi=(sites[i]['y']+dy)//ly;wj=(sites[j]['y']+dy)//ly
  key=bkey(sm[i],sm[j],b['wy']+wj-wi)
  assert key in bl
  bm[k]=bl[key]
 assert len(set(sm.values()))==len(ix) and len(set(bm.values()))==len(ib)
 return sm,bm

def calc(a,b,ix,ib,sm,bm):
 def stat(v):return {'rms':math.sqrt(math.fsum(x*x for x in v)/len(v)),'max_abs':max(map(abs,v))}
 return {'sz':stat([b['sz_profile'][sm[i]]-a['sz_profile'][i] for i in ix]),'bond':stat([b['bond_energy'][bm[k]]-a['bond_energy'][k] for k in ib])}
answers={}
for left,right in itertools.combinations(('random','hourglass','windmill'),2):
 a,b=rows[left],rows[right];entries=[]
 for dx in (-1,0,1):
  ix,ib=support(dx)
  baseline=calc(a,b,ix,ib,{i:i for i in ix},{k:k for k in ib})
  yonly=[{'dy':dy,**calc(a,b,ix,ib,*mapping(0,dy,ix,ib))} for dy in range(3)]
  besty=min(yonly,key=lambda r:r['bond']['rms'])
  for dy in range(3):
   value=calc(a,b,ix,ib,*mapping(dx,dy,ix,ib))
   entries.append({'dx':dx,'dy':dy,'source_x':sorted({sites[i]['x'] for i in ix}),'target_x':sorted({sites[i]['x']+dx for i in ix}),
    'site_count':len(ix),'bond_count':len(ib),**value,'same_source_support_dx0_dy0':baseline,'best_y_only_same_source_support':besty,
    'bond_ratio_to_same_support_raw':value['bond']['rms']/baseline['bond']['rms'] if baseline['bond']['rms'] else None,
    'bond_ratio_to_best_y_only_same_support':value['bond']['rms']/besty['bond']['rms'] if besty['bond']['rms'] else None})
 answers[f'{left}__{right}']={'all_shifts':entries,'best_bond_rms_posthoc':min(entries,key=lambda e:e['bond']['rms'])}
# Recheck the immutable completed input and selected states before publication.
now=tomllib.loads(p.read_text());newrows={r['branch']:r for r in now['batches'] if r['stage']==args.stage}
assert rows==newrows and now['configuration']==g and p.read_bytes()==raw
out={'scope':args.stage+'_posthoc_axial_domain_diagnostic','analysis_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'worker_status_at_read':d['status'],'validation_sha256_at_read':hashlib.sha256(raw).hexdigest(),
 'selected_rows_and_input_bytes_unchanged_at_end':True,'all_selected_checkpoint_pins_verified':True,'axial_translation_is_whole_system_symmetry':False,
 'definition':'compare left(x,y,s) with right(x+dx,y+dy,s); source and target both inside x=1..4; internal bonds only',
 'support_warning':'dx0 has36sites/66bonds; dx nonzero has27sites/48bonds. Cross-support minimum is descriptive; same-support dx0 and best-y-only baselines are included.',
 'pairs':answers}
args.output.parent.mkdir(parents=True,exist_ok=True)
with args.output.open('x') as f: json.dump(out,f,indent=2,allow_nan=False);f.write('\n')
for pair,a in answers.items():
 print(pair)
 for dx in (-1,1):
  e=min((e for e in a['all_shifts'] if e['dx']==dx),key=lambda e:e['bond']['rms'])
  print('dx',dx,'best dy',e['dy'],'bond RMS',e['bond']['rms'],'same-support y-only ratio',e['bond_ratio_to_best_y_only_same_support'])
 print('best',a['best_bond_rms_posthoc']['dx'],a['best_bond_rms_posthoc']['dy'])
