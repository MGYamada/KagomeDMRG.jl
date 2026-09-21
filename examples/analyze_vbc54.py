#!/usr/bin/env python3
"""Audit and compare the completed three-branch N54 preparation experiment.

Standard-library analysis; matplotlib is imported only when rendering figures.
No DMRG, MPS loading, phase assignment, or energy uncertainty extrapolation.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
import tomllib

ROOT = Path(__file__).resolve().parent.parent
BRANCHES = ('random', 'hourglass', 'windmill')
STAGES = ('prepare', 'reduce', 'release', 'relax')
LIMITS = dict(energy_per_site=1e-6, sz_profile=1e-4, bond_profile=1e-4,
              variance_per_site=1e-5, truncation=1e-6)


def require(test, message):
    if not test:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def label(path):
    path = Path(path).resolve()
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def read_toml(path):
    return tomllib.loads(Path(path).read_text())


def rootpath(path):
    p = Path(path)
    return p.resolve() if p.is_absolute() else (ROOT / p).resolve()


def finite(values):
    return all(isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x) for x in values)


def stats(values):
    require(bool(values) and finite(values), 'empty or nonfinite observable')
    mean = math.fsum(values) / len(values)
    return dict(count=len(values), mean=mean,
                rms=math.sqrt(math.fsum(x*x for x in values)/len(values)),
                centered_rms=math.sqrt(math.fsum((x-mean)**2 for x in values)/len(values)),
                max_abs=max(map(abs, values)), minimum=min(values), maximum=max(values))


def difference(left, right, indices):
    return stats([right[i]-left[i] for i in indices])


def bond_key(i, j, wy):
    return (i, j, wy) if i < j else (j, i, -wy)


def shift_permutations(geometry, shift):
    """Pull back right profiles: at site (x,y,s), sample right(x,y+shift,s).

    Endpoint wrap integers n_i,n_j change winding to wy+n_j-n_i.
    Orienting the translated bond backwards also reverses its winding.
    """
    ly = geometry['Ly']
    sites, bonds = geometry['sites'], geometry['bonds']
    lookup = {(s['x'], s['y'], s['sublattice']): i for i, s in enumerate(sites)}
    mapping = [lookup[(s['x'], (s['y']+shift) % ly, s['sublattice'])] for s in sites]
    images = [(s['y']+shift)//ly for s in sites]
    bond_lookup = {bond_key(b['i'], b['j'], b['wy']): k for k, b in enumerate(bonds)}
    require(len(bond_lookup) == len(bonds), 'duplicate physical bonds')
    mapped = []
    for b in bonds:
        i, j = b['i']-1, b['j']-1
        key = bond_key(mapping[i]+1, mapping[j]+1, b['wy']+images[j]-images[i])
        require(key in bond_lookup, 'translation does not preserve bond geometry/winding')
        mapped.append(bond_lookup[key])
    require(sorted(mapping) == list(range(len(sites))) and sorted(mapped) == list(range(len(bonds))),
            'translation is not bijective')
    return mapping, mapped


def bond_descriptor(geometry, bond):
    sites = geometry['sites']
    a, b = sites[bond['i']-1], sites[bond['j']-1]
    wy = bond['wy']
    if a['sublattice'] > b['sublattice']:
        a, b, wy = b, a, -wy
    dx, dy = b['x']-a['x'], b['y']+geometry['Ly']*wy-a['y']
    kind = f"{a['sublattice']}{b['sublattice']}({dx},{dy})"
    return kind, a['x'], a['y']


def window(geometry, edge):
    keep = [i for i, s in enumerate(geometry['sites']) if edge <= s['x'] < geometry['Lx']-edge]
    selected = set(keep)
    bonds = [i for i, b in enumerate(geometry['bonds']) if b['i']-1 in selected and b['j']-1 in selected]
    return keep, bonds


def projection(geometry, row, weights, indices):
    """Remove each of the six oriented NN-kind means from both arrays."""
    groups = {}
    for i in indices:
        groups.setdefault(bond_descriptor(geometry, geometry['bonds'][i])[0], []).append(i)
    values, template = [], []
    for group in groups.values():
        m = math.fsum(row['bond_energy'][i] for i in group)/len(group)
        w = math.fsum(weights[i] for i in group)/len(group)
        values.extend(row['bond_energy'][i]-m for i in group)
        template.extend(weights[i]-w for i in group)
    dot = math.fsum(x*y for x, y in zip(values, template))
    vv, ww = math.fsum(x*x for x in values), math.fsum(x*x for x in template)
    return dict(dot_per_bond=dot/len(indices), amplitude=dot/ww if ww > 0 else None,
                cosine=dot/math.sqrt(vv*ww) if vv*ww > 0 else None,
                bond_kind_mean_removed=True, template_squared_norm=ww,
                interpretation='geometric template overlap, not a phase or wavefunction fidelity')


def profile_summary(geometry, row, templates, edge):
    sites, bonds = window(geometry, edge)
    sub = {s: stats([row['sz_profile'][i] for i in sites if geometry['sites'][i]['sublattice'] == s])
           for s in 'ABC'}
    kinds = sorted({bond_descriptor(geometry, geometry['bonds'][i])[0] for i in bonds})
    return dict(edge_columns_excluded=edge, included_columns=list(range(edge, geometry['Lx']-edge)),
                bond_selection='both endpoints in retained columns',
                sz=stats([row['sz_profile'][i] for i in sites]),
                bond=stats([row['bond_energy'][i] for i in bonds]), sz_by_sublattice=sub,
                bond_by_kind={k: stats([row['bond_energy'][i] for i in bonds
                                       if bond_descriptor(geometry, geometry['bonds'][i])[0] == k]) for k in kinds},
                template_projection={name: projection(geometry, row, t['weights'], bonds)
                                     for name, t in templates.items()})


def pair_comparison(geometry, left, right, edge):
    sites, bonds = window(geometry, edge)
    shifts = []
    for y in range(geometry['Ly']):
        smap, bmap = shift_permutations(geometry, y)
        shifts.append(dict(y_shift=y,
            bond=difference(left['bond_energy'], [right['bond_energy'][i] for i in bmap], bonds),
            sz=difference(left['sz_profile'], [right['sz_profile'][i] for i in smap], sites)))
    best = min(shifts, key=lambda r: (r['bond']['rms'], r['y_shift']))
    return dict(difference_convention='translated right minus left', raw=shifts[0],
                all_circumference_shifts=shifts, minimum_bond_rms_shift=best,
                sz_alignment='Sz evaluated at the same shift minimizing bond RMS', axial_shift='not_tested')


def correlation_lines(geometry, row):
    sites = geometry['sites']
    ids = {(s['x'], s['y'], s['sublattice']): i for i, s in enumerate(sites)}
    anchor = ids[(2, 0, 'A')]
    result = dict(anchor=dict(index=anchor+1, x=2, y=0, sublattice='A'), fixed_Q_Splus_mean=0.0)
    for direction, coords in [('axis', [(x, 0, 'A') for x in range(6)]),
                              ('circumference', [(2, y, 'A') for y in range(3)])]:
        rows = []
        for coord in coords:
            j = ids[coord]
            rows.append(dict(index=j+1, x=coord[0], y=coord[1], sublattice='A',
                cell_displacement=coord[0]-2 if direction == 'axis' else coord[1],
                connected_zz=row['correlations']['zz_real'][anchor][j]-row['sz_profile'][anchor]*row['sz_profile'][j],
                pm_real=row['correlations']['pm_real'][anchor][j],
                pm_imag=row['correlations']['pm_imag'][anchor][j]))
        result[direction] = rows
    result['interpretation'] = 'single-anchor same-sublattice paths; no correlation-length or gap fit'
    return result


def audit_observables(geometry, row):
    n = geometry['N']
    sz, bonds, corr = row['sz_profile'], row['bond_energy'], row['correlations']
    require(len(sz) == n and len(bonds) == len(geometry['bonds']) and finite(sz+bonds), 'profile shape/finite failure')
    for key in ('zz_real', 'zz_imag', 'pm_real', 'pm_imag'):
        require(len(corr[key]) == n and all(len(r) == n and finite(r) for r in corr[key]), 'correlation matrix malformed')
    errors = dict(charge=abs(math.fsum(sz)-3), bond_sum=abs(math.fsum(bonds)-row['energy']),
                  variance=abs(row['HdaggerH_real']-row['energy']**2-row['variance']))
    errors['bond_from_correlations'] = max(abs(bonds[k]-corr['zz_real'][b['i']-1][b['j']-1]-corr['pm_real'][b['i']-1][b['j']-1])
                                           for k,b in enumerate(geometry['bonds']))
    errors['fixed_charge_zz'] = max(abs(math.fsum(corr['zz_real'][i])-3*sz[i]) for i in range(n))
    errors['pm_hermiticity'] = max(abs(complex(corr['pm_real'][i][j]-corr['pm_real'][j][i],
                                                           corr['pm_imag'][i][j]+corr['pm_imag'][j][i])) for i in range(n) for j in range(n))
    require(max(errors.values()) <= 1e-9, f'independent sum-rule check failed: {errors}')
    require(row['variance'] >= -row['variance_roundoff_scale'] and abs(row['HdaggerH_imag']) <= row['variance_roundoff_scale'],
            'invalid variance/complex second moment')
    return errors


def load_audited(path):
    path = Path(path).resolve()
    record = read_toml(path)
    require(record.get('status') == 'completed_preparation_comparison_accuracy_separate', 'worker is incomplete')
    require(record.get('schema_version') == 1 and record.get('N') == 54 and record.get('Q') == 6 and record.get('theta') == 0,
            'wrong schema/model sector')
    require(record.get('sources_unchanged') is True and record.get('backend_unchanged') is True, 'source finalization missing')
    run = path.parent
    execution_path = run/'execution.toml'
    execution = read_toml(execution_path)
    require(execution.get('status') == 'exited' and execution.get('worker_exit_confirmed') is True
            and execution.get('worker_exit_code') == 0, 'worker execution did not complete successfully')
    require(execution['julia_threads'] == execution['blas_threads'] == record['runtime']['julia_threads'] == record['runtime']['blas_threads'] == 1,
            'thread allocation mismatch')
    require(0 <= execution['wall_elapsed_seconds'] <= execution['wall_limit_seconds'] == record['wall_limit_seconds'] == 2700,
            'wall budget mismatch')
    archive = run/'analysis-sources'
    pins = {label(path): digest(path), label(execution_path): digest(execution_path)}
    cfgpath = archive/'config.toml'
    require(digest(cfgpath) == record['config_sha256'] and read_toml(cfgpath) == record['config'], 'config archive mismatch')
    pins[label(cfgpath)] = digest(cfgpath)
    require(record['config']['precision'] == LIMITS and record['config']['branches'] == list(BRANCHES), 'protocol mismatch')
    for name, expected in record['analysis_source_sha256'].items():
        source = archive/name
        require(digest(source) == expected, f'archived source mismatch: {name}')
        pins[label(source)] = expected
    geometry = record['configuration']
    require((geometry['Lx'],geometry['Ly'],geometry['N'],geometry['Q']) == (6,3,54,6)
            and geometry['gauge'] == 'seam' and geometry['ordering'] == 'x_then_y_then_A_B_C'
            and all(h == 0 for h in geometry['hz']) and all(b['Jxy'] == b['Jz'] == 1 for b in geometry['bonds']), 'wrong final Hamiltonian')
    require(len(record['batches']) == 12, 'missing preparation stages')
    require([(r['stage'], r['branch']) for r in record['batches']] == list(itertools.product(STAGES, BRANCHES)),
            'stage ordering or branch set mismatch')
    require(len(geometry['sites']) == 54 and len(geometry['bonds']) == 102, 'wrong geometry sizes')
    selected = {}
    checks = {}
    common_sites = None
    for row in record['batches']:
        cp = rootpath(row['checkpoint'])
        require(cp.is_relative_to(run), 'checkpoint is not bound to the run')
        require(row.get('checkpoint_verified') is True and row['checkpoint_overlap_error'] <= 1e-12, 'checkpoint reload failed')
        actual = {name: digest(cp/name) for name in ('metadata.toml', 'state.jls', 'checksums.toml')}
        require(actual == row['checkpoint_sha256'], 'checkpoint pin mismatch')
        checksum = read_toml(cp/'checksums.toml')
        for name in ('metadata.toml', 'state.jls'):
            require(checksum['files'][name] == dict(sha256=actual[name], bytes=(cp/name).stat().st_size), 'checkpoint checksum mismatch')
        pins.update({label(cp/name): sha for name,sha in actual.items()})
        meta = read_toml(cp/'metadata.toml')
        require(meta['status'] == 'trial' and meta['theta_path'] == [0.0] and meta['configuration'] == row['model_configuration'],
                'checkpoint model/status mismatch')
        require(meta['provenance']['source_sha256'] == record['code']['source_sha256'] and meta['settings'] == row['settings'],
                'checkpoint source/settings mismatch')
        require(meta['provenance']['environment_sha256'] == record['code']['environment_sha256']
                and meta['provenance']['manifest'] == record['code']['manifest'], 'checkpoint environment mismatch')
        require(meta['runtime'] == {k:v for k,v in record['runtime'].items() if k not in ('julia_threads','blas_threads')}, 'runtime mismatch')
        state = meta['state']
        if common_sites is None:
            common_sites = state['site_indices']
        require(state['site_indices'] == common_sites, 'branches or stages changed physical site indices')
        stage = record['config']['stages'][STAGES.index(row['stage'])]
        expected_lambda = 0.0 if row['branch'] == 'random' else stage['lambda']
        require(row['lambda'] == expected_lambda and row['settings']['nsweeps'] == 2
                and row['settings']['maxdim'] == stage['maxdim'], 'actual preparation schedule mismatch')
        model = row['model_configuration']
        require({k:v for k,v in model.items() if k != 'bonds'} ==
                {k:v for k,v in geometry.items() if k != 'bonds'}, 'preparation changed non-bond model data')
        weights = [0.0]*102 if row['branch'] == 'random' else record['templates'][row['branch']]['weights']
        require(len(model['bonds']) == len(geometry['bonds']) == len(weights), 'prepared bond count mismatch')
        for actual_bond, base_bond, weight in zip(model['bonds'],geometry['bonds'],weights):
            require({k:v for k,v in actual_bond.items() if k not in ('Jxy','Jz')} ==
                    {k:v for k,v in base_bond.items() if k not in ('Jxy','Jz')}, 'preparation changed bond geometry')
            expected_J = 1+expected_lambda*weight
            require(abs(actual_bond['Jxy']-expected_J) <= 1e-14 and actual_bond['Jxy'] == actual_bond['Jz'],
                    'prepared coupling does not implement the saved template')
        require(state['Q'] == 6 and state['theta'] == 0 and state['energy'] == row['energy'] and state['sz'] == row['sz_profile']
                and state['sweep_energies'] == row['sweep_energies'] and state['max_truncation_errors'] == row['measured_truncation_errors'], 'checkpoint observable mismatch')
        require(abs(abs(state['local_energy']-state['energy'])-row['last_optimizer_energy_error']) <= 1e-12, 'optimizer error mismatch')
        require(len(row['sweep_energies']) == 2 and len(row['measured_truncation_errors']) == 2
                and all(0 <= x <= 1 for x in row['measured_truncation_errors'])
                and abs(abs(row['sweep_energies'][-1]-row['sweep_energies'][-2])-row['last_sweep_energy_change']) <= 1e-12,
                'sweep count or last-sweep diagnostic mismatch')
        if row['stage'] in ('release', 'relax'):
            require(row.get('integrity_passed') is True and row.get('active_phase') == 'completed'
                    and row['status'] == 'diagnostics_passed_accuracy_separate' and row['lambda'] == 0
                    and row['final_nn_model_verified'] is True and row['model_configuration'] == geometry,
                    'incomplete or pinned final row')
            checks[f"{row['branch']}_{row['stage']}"] = audit_observables(geometry, row)
            selected[(row['branch'],row['stage'])] = row
    for name in ('hourglass','windmill'):
        require(len(record['templates'][name]['weights']) == 102 and finite(record['templates'][name]['weights']), 'template weights malformed')
    for y in range(3):
        shift_permutations(geometry, y)
    return record, selected, pins, checks, execution


def analyze(record, selected):
    geometry = record['configuration']
    result = dict(schema_version=1, scope='N54_Q6_finite_chi_preparation_sensitivity',
        model='nearest-neighbor isotropic J=1; theta=0; final preparation lambda=0',
        geometry=geometry, templates=record['templates'],
        variance_interpretation='Hamiltonian dispersion, not a ground-energy error bar',
        phase_identification='not_attempted', axial_alignment='not_tested', branches={}, pairs={})
    for branch in BRANCHES:
        before, final = selected[(branch,'release')], selected[(branch,'relax')]
        require(before['cumulative_sweeps'] == 6 and final['cumulative_sweeps'] == 8, 'unexpected sweep count')
        values = dict(energy_per_site=max(abs(final['energy']-before['energy']),final['last_sweep_energy_change'],final['last_optimizer_energy_error'])/54,
            sz_profile=max(abs(a-b) for a,b in zip(before['sz_profile'],final['sz_profile'])),
            bond_profile=max(abs(a-b) for a,b in zip(before['bond_energy'],final['bond_energy'])),
            variance_per_site=abs(final['variance'])/54, truncation=final['measured_truncation_errors'][-1])
        precision = {k:dict(value=v, limit=LIMITS[k], passed=v <= LIMITS[k]) for k,v in values.items()}
        require(all(abs(final['precision_values'][k]-v) <= 1e-12 and final['precision_passed'][k] == precision[k]['passed'] for k,v in values.items())
                and final['all_precision_conditions_passed'] == all(v['passed'] for v in precision.values()), 'worker precision arithmetic mismatch')
        summaries = {}
        for stage, row in [('release6',before),('relax8',final)]:
            summaries[stage] = dict(energy=row['energy'], energy_per_site=row['energy']/54,
                variance=row['variance'], variance_per_site=row['variance']/54,
                measured_truncation_error=row['measured_truncation_errors'][-1], maxlinkdim=row['maxlinkdim'],
                sz_profile=row['sz_profile'], bond_energy=row['bond_energy'],
                windows={str(edge):profile_summary(geometry,row,record['templates'],edge) for edge in (1,2)},
                correlations=correlation_lines(geometry,row))
        summaries['release6_to_relax8'] = dict(energy_change=final['energy']-before['energy'],
            windows={str(edge):dict(sz=difference(before['sz_profile'],final['sz_profile'],window(geometry,edge)[0]),
                                  bond=difference(before['bond_energy'],final['bond_energy'],window(geometry,edge)[1])) for edge in (1,2)},
            precision=precision, all_precision_conditions_passed=all(v['passed'] for v in precision.values()))
        result['branches'][branch] = summaries
    for a,b in itertools.combinations(BRANCHES,2):
        result['pairs'][f'{a}__{b}'] = {stage:dict(energy_difference=selected[(b,stage)]['energy']-selected[(a,stage)]['energy'],
            windows={str(edge):pair_comparison(geometry,selected[(a,stage)],selected[(b,stage)],edge) for edge in (1,2)})
            for stage in ('release','relax')}
    result['claim_status'] = 'all_selected_precision_checks_passed_but_chi_size_seed_dependence_unestablished' if all(
        result['branches'][b]['release6_to_relax8']['all_precision_conditions_passed'] for b in BRANCHES) else 'unconverged_trials_no_energy_ranking_claim'
    return result


def render(result, prefix):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    geometry = result['geometry']
    finals = [result['branches'][b]['relax8'] for b in BRANCHES]
    szlim = max(abs(x) for r in finals for x in r['sz_profile']) or 1e-12
    blo = min(x for r in finals for x in r['bond_energy'])
    bhi = max(x for r in finals for x in r['bond_energy'])
    descriptors = [bond_descriptor(geometry,b) for b in geometry['bonds']]
    kinds = sorted({d[0] for d in descriptors})
    fig, axes = plt.subplots(4,3,figsize=(15,14),layout='constrained',sharey='row')
    for col,(branch,row) in enumerate(zip(BRANCHES,finals)):
        szgrid = [[math.nan]*6 for _ in range(9)]
        for i,s in enumerate(geometry['sites']):
            szgrid[3*s['y']+'ABC'.index(s['sublattice'])][s['x']] = row['sz_profile'][i]
        ax=axes[0,col];im=ax.imshow(szgrid,aspect='auto',origin='lower',cmap='RdBu_r',vmin=-szlim,vmax=szlim)
        ax.set_title(f"{branch} | E/J={row['energy']:.8f}",fontsize=11)
        ax.set_yticks(range(9),[f'{s}, y={y}' for y in range(3) for s in 'ABC']);ax.set_xlabel('column x')
        grid = [[math.nan]*6 for _ in range(3*len(kinds))]
        for value,(kind,x,y) in zip(row['bond_energy'],descriptors):grid[3*kinds.index(kind)+y][x]=value
        ax=axes[1,col];bm=ax.imshow(grid,aspect='auto',origin='lower',cmap='viridis',vmin=blo,vmax=bhi)
        ax.set_yticks([3*i+1 for i in range(len(kinds))],kinds,fontsize=8);ax.set_xlabel('bond anchor column x')
        for r in (0,1):
            axes[r,col].axvline(.5,color='white',ls='--',lw=.8);axes[r,col].axvline(4.5,color='white',ls='--',lw=.8)
            axes[r,col].axvspan(1.5,3.5,facecolor='none',edgecolor='#f4a340',lw=1.3)
        for r,direction in ((2,'axis'),(3,'circumference')):
            ax=axes[r,col];line=[v for v in row['correlations'][direction] if v['cell_displacement'] != 0]
            x=[v['cell_displacement'] for v in line]
            for key,name,style in [('connected_zz','connected SzSz','o-'),('pm_real','Re S+S-','s-'),('pm_imag','Im S+S-','^-')]:
                ax.plot(x,[v[key] for v in line],style,label=name,markersize=4)
            ax.axhline(0,color='0.5',lw=.6);ax.set_xlabel('axial cell displacement' if direction=='axis' else 'circumference cell displacement')
            ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.colorbar(im,ax=axes[0,:],label='physical Sz (common scale)',shrink=.85)
    fig.colorbar(bm,ax=axes[1,:],label='bond energy / J (common scale)',shrink=.85)
    title='UNCONVERGED TRIALS' if result['claim_status'].startswith('unconverged') else 'FINITE-CHI TRIALS'
    fig.suptitle(f'N54, Q6, unpinned nearest-neighbor model | {title}\n'
                 'Total 8 sweeps; NN relaxation: random 8, prepared 4 | wrap 3a2\n'
                 'Dashed: exclude 1 edge column; orange: anchor columns x=2,3\n'
                 'Off-site correlations: A(2,0) to A sites; no gap or phase inference',fontsize=13)
    axes[1,0].set_ylabel('NN kind; each triplet: y=0,1,2 bottom to top',fontsize=8)
    prefix=Path(prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    for extension in ('png','svg'):
        target=prefix.with_suffix('.'+extension)
        require(not target.exists(),f'refusing existing figure: {target}')
        fig.savefig(target,dpi=180)
    plt.close(fig)
    return dict(matplotlib_version=matplotlib.__version__,
                figure_sha256={label(prefix.with_suffix('.'+ext)):digest(prefix.with_suffix('.'+ext))
                               for ext in ('png','svg')})


def shift_self_test():
    # Independent infinite-lattice NN templates; includes boundary-crossing bonds.
    sites=[dict(index=1+9*x+3*y+s,x=x,y=y,sublattice='ABC'[s]) for x in range(6) for y in range(3) for s in range(3)]
    index={(s['x'],s['y'],s['sublattice']):s['index'] for s in sites}
    bonds=[]
    for x in range(6):
        for y in range(3):
            for a,b,dx,dy in [('A','B',0,0),('A','C',0,0),('B','C',0,0),('A','B',-1,0),('A','C',0,-1),('B','C',1,-1)]:
                if 0<=x+dx<6:
                    i,j=index[(x,y,a)],index[(x+dx,(y+dy)%3,b)]
                    bonds.append(dict(i=i,j=j,wy=(y+dy)//3))
    geometry=dict(Lx=6,Ly=3,N=54,sites=sites,bonds=bonds)
    require(len(bonds)==102,'fixture bond count')
    sitevalues=[float(i) for i in range(54)];bondvalues=[float(i*i+3*i) for i in range(102)]
    sm,bm=shift_permutations(geometry,1)
    inverse_s=[0.0]*54;inverse_b=[0.0]*102
    for i,j in enumerate(sm):inverse_s[j]=sitevalues[i]
    for i,j in enumerate(bm):inverse_b[j]=bondvalues[i]
    answer=pair_comparison(geometry,dict(sz_profile=sitevalues,bond_energy=bondvalues),dict(sz_profile=inverse_s,bond_energy=inverse_b),1)
    require(answer['minimum_bond_rms_shift']['y_shift']==1 and answer['minimum_bond_rms_shift']['bond']['max_abs']==0
            and answer['minimum_bond_rms_shift']['sz']['max_abs']==0,'synthetic shifted pattern not recovered')
    require(all(bm[bm[bm[i]]]==i for i in range(102)),'three translations do not close')
    # Reversing every bond must leave energy translation unchanged.
    reversed_geometry={**geometry,'bonds':[dict(i=b['j'],j=b['i'],wy=-b['wy']) for b in bonds]}
    require(shift_permutations(reversed_geometry,1)==(sm,bm),'orientation reversal changes scalar-energy mapping')
    bad={**geometry,'bonds':[dict(b) for b in bonds]}
    bad['bonds'][0]['wy'] += 1
    try:
        shift_permutations(bad,1)
    except ValueError:
        pass
    else:
        raise ValueError('incorrect winding was not rejected')
    return dict(status='passed',fixture='independent_NN_templates_102_bonds',checks=['shift1 recovery','three-cycle closure','reversed winding/orientation','wrong winding rejected'],not_a_DMRG_result=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('validation',type=Path,nargs='?')
    parser.add_argument('--output',type=Path,default=ROOT/'docs/research/data/p4_vbc54_comparison.json')
    parser.add_argument('--figure-prefix',type=Path,default=ROOT/'docs/research/figures/p4_vbc54_profiles')
    parser.add_argument('--no-plot',action='store_true')
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args()
    if args.self_test:
        print(json.dumps(shift_self_test(),indent=2))
        if args.validation is None:return
    if args.validation is None:parser.error('validation is required unless --self-test')
    require(not args.output.exists(),f'refusing existing output: {args.output}')
    if not args.no_plot:
        for ext in ('png','svg'):require(not args.figure_prefix.with_suffix('.'+ext).exists(),'refusing existing figure')
    record,selected,pins,checks,execution=load_audited(args.validation)
    result=analyze(record,selected)
    result['provenance']=dict(input_sha256=pins,analysis_path=label(__file__),analysis_sha256=digest(__file__),
        python_version=sys.version.split()[0],
        checks=checks,shift_fixture=shift_self_test(),sources_unchanged=record['sources_unchanged'],backend_unchanged=record['backend_unchanged'],
        checkpoint_payloads_hash_verified=True,mps_recontracted=False,
        execution={k:execution[k] for k in ('status','wall_elapsed_seconds','wall_limit_seconds','worker_exit_code','worker_exit_confirmed','julia_threads','blas_threads')})
    if not args.no_plot:result['artifacts']=render(result,args.figure_prefix)
    # Recheck all input bytes after computation/plotting before publishing analysis.
    require(all(digest(rootpath(p))==h for p,h in pins.items()),'input changed during analysis')
    result['provenance']['inputs_unchanged_after_analysis']=True
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:json.dump(result,out,ensure_ascii=False,indent=2,allow_nan=False);out.write('\n')
    print(json.dumps(dict(output=label(args.output),status=result['claim_status'],branches={b:result['branches'][b]['release6_to_relax8']['precision'] for b in BRANCHES}),indent=2))


if __name__=='__main__':
    try:main()
    except (ValueError,KeyError,OSError) as error:
        raise SystemExit(str(error)) from error
