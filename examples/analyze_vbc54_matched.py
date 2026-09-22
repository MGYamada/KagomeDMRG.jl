#!/usr/bin/env python3
"""Audit N54 matched-parent children and compare structure after equal sweeps.

No DMRG or MPS deserialization. Every correlation comparison translates both
site indices, excludes on-site terms, and uses the bond-selected alignment.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import analyze_vbc54 as base

ROOT = base.ROOT
BRANCHES = ('random', 'windmill')
CHIS = (128, 256)
EDGES = (1, 2)
FILES = ('metadata.toml', 'state.jls', 'checksums.toml')
require, digest, label = base.require, base.digest, base.label
read_toml, rootpath = base.read_toml, base.rootpath


def audit_record(path, status, budget, pins):
    path = Path(path).resolve()
    record = read_toml(path)
    require(record.get('status') == status, 'worker is incomplete')
    require((record['schema_version'], record['N'], record['Q'], record['theta']) == (1, 54, 6, 0),
            'wrong schema/model sector')
    require(record.get('sources_unchanged') is True and record.get('backend_unchanged') is True,
            'source finalization missing')
    ep = path.parent / 'execution.toml'
    execution = read_toml(ep)
    require(execution.get('status') == 'exited' and execution.get('worker_exit_confirmed') is True
            and execution.get('worker_exit_code') == 0, 'worker execution incomplete')
    require(execution['julia_threads'] == execution['blas_threads'] ==
            record['runtime']['julia_threads'] == record['runtime']['blas_threads'] == 1,
            'thread allocation mismatch')
    require(0 <= execution['wall_elapsed_seconds'] <= execution['wall_limit_seconds'] ==
            record['wall_limit_seconds'] == record['config']['wall_seconds'] == budget, 'wall budget mismatch')
    archive = path.parent / 'analysis-sources'
    cfg = archive / 'config.toml'
    require(digest(cfg) == record['config_sha256'] and read_toml(cfg) == record['config'],
            'archived configuration mismatch')
    pins.update({label(path): digest(path), label(ep): digest(ep), label(cfg): digest(cfg)})
    for name, sha in record['analysis_source_sha256'].items():
        source = (archive / name).resolve()
        require(source.is_relative_to(archive), 'archived source path escaped archive')
        require(digest(source) == sha, f'archived source mismatch: {name}')
        pins[label(source)] = sha
    g = record['configuration']
    require((g['Lx'], g['Ly'], g['N'], g['Q']) == (6, 3, 54, 6)
            and g['gauge'] == 'seam' and g['ordering'] == 'x_then_y_then_A_B_C'
            and len(g['sites']) == 54 and len(g['bonds']) == 102
            and all(h == 0 for h in g['hz']) and all(b['Jxy'] == b['Jz'] == 1 for b in g['bonds']),
            'wrong final Hamiltonian')
    require(record['config']['precision'] == base.LIMITS, 'changed precision criteria')
    for shift in range(3):
        base.shift_permutations(g, shift)
    return record, execution


def audit_checkpoint(row, record, run, pins):
    cp = rootpath(row['checkpoint'])
    require(cp.is_relative_to(run), 'checkpoint is not bound to its run')
    require(row.get('checkpoint_verified') is True and row['checkpoint_overlap_error'] <= 1e-12,
            'checkpoint reload failed')
    actual = {name: digest(cp / name) for name in FILES}
    require(actual == row['checkpoint_sha256'], 'checkpoint pin mismatch')
    checksums = read_toml(cp / 'checksums.toml')
    for name in FILES[:2]:
        require(checksums['files'][name] == dict(sha256=actual[name], bytes=(cp / name).stat().st_size),
                'checkpoint checksum mismatch')
    pins.update({label(cp / name): sha for name, sha in actual.items()})
    meta = read_toml(cp / 'metadata.toml')
    require(meta['status'] == 'trial' and meta['theta_path'] == [0.0]
            and meta['configuration'] == row['model_configuration'] == record['configuration'],
            'checkpoint model/status mismatch')
    require(meta['settings'] == row['settings'], 'checkpoint settings mismatch')
    for key in ('source_sha256', 'environment_sha256', 'manifest'):
        require(meta['provenance'][key] == record['code'][key], f'checkpoint {key} mismatch')
    require(meta['runtime'] == {k: v for k, v in record['runtime'].items()
                               if k not in ('julia_threads', 'blas_threads')}, 'runtime mismatch')
    state = meta['state']
    require(state['Q'] == 6 and state['theta'] == 0 and state['energy'] == row['energy']
            and state['sz'] == row['sz_profile'] and state['sweep_energies'] == row['sweep_energies']
            and state['max_truncation_errors'] == row['measured_truncation_errors'],
            'checkpoint observable mismatch')
    require(len(row['sweep_energies']) == len(row['measured_truncation_errors']) == 2
            and base.finite(row['sweep_energies'] + row['measured_truncation_errors'])
            and all(0 <= e <= 1 for e in row['measured_truncation_errors']), 'bad sweep diagnostics')
    require(abs(abs(state['local_energy'] - row['energy']) - row['last_optimizer_energy_error']) <= 1e-12
            and abs(abs(row['sweep_energies'][1] - row['sweep_energies'][0]) -
                    row['last_sweep_energy_change']) <= 1e-12, 'last-sweep diagnostic mismatch')
    require(row.get('integrity_passed') is True and row['active_phase'] == 'completed'
            and row['status'] == 'diagnostics_passed_accuracy_separate' and row['lambda'] == 0
            and row['final_nn_model_verified'] is True, 'incomplete or pinned row')
    require(all(math.isfinite(value) and 0 <= value <= record['integrity_limits'][key]
                for key, value in row['integrity_errors'].items()), 'recorded integrity check failed')
    require(abs(row['energy_per_site'] - row['energy'] / 54) <= 1e-12
            and abs(row['variance_per_site'] - row['variance'] / 54) <= 1e-12, 'per-site arithmetic mismatch')
    return meta, base.audit_observables(record['configuration'], row)


def precision(parent, child):
    values = dict(energy_per_site=max(abs(child['energy'] - parent['energy']),
                  child['last_sweep_energy_change'], child['last_optimizer_energy_error']) / 54,
        sz_profile=max(abs(a-b) for a, b in zip(parent['sz_profile'], child['sz_profile'])),
        bond_profile=max(abs(a-b) for a, b in zip(parent['bond_energy'], child['bond_energy'])),
        variance_per_site=abs(child['variance']) / 54,
        truncation=child['measured_truncation_errors'][-1])
    result = {k: dict(value=v, limit=base.LIMITS[k], passed=v <= base.LIMITS[k]) for k, v in values.items()}
    require(set(child['precision_values']) == set(values) == set(child['precision_passed'])
            and all(abs(child['precision_values'][k]-v) <= 1e-12
                    and child['precision_passed'][k] == result[k]['passed'] for k, v in values.items())
            and child['all_precision_conditions_passed'] == all(r['passed'] for r in result.values()),
            'worker precision arithmetic mismatch')
    return result


def load_audited(path):
    path = Path(path).resolve()
    pins, checks = {}, {}
    record, execution = audit_record(path, 'completed_matched_parent_comparison_accuracy_separate', 3600, pins)
    cfg = record['config']
    require(record.get('parents_unchanged') is True and cfg['branches'] == list(BRANCHES)
            and cfg['maxdims'] == list(CHIS) and cfg['additional_sweeps'] == 2, 'changed matched protocol')
    require([(r['chi'], r['branch']) for r in record['batches']] == list(itertools.product(CHIS, BRANCHES)),
            'missing/reordered/extra matched children')
    require(len(cfg['parents']) == 2 and [p['branch'] for p in cfg['parents']] == list(BRANCHES)
            and set(record['parents']) == set(BRANCHES), 'wrong parent set')
    parents, metas, templates, cache = {}, {}, None, {}
    for spec in cfg['parents']:
        branch = spec['branch']
        pp = rootpath(spec['record'])
        require(digest(pp) == spec['record_sha256'], 'parent record hash mismatch')
        if pp not in cache:
            cache[pp] = audit_record(pp, 'completed_preparation_comparison_accuracy_separate', 2700, pins)[0]
        pr = cache[pp]
        require(pr['configuration'] == record['configuration'], 'parent geometry/model mismatch')
        for key in ('source_sha256', 'environment_sha256', 'manifest'):
            require(pr['code'][key] == record['code'][key], f'parent/child {key} mismatch')
        require(pr['runtime'] == record['runtime'], 'parent/child runtime mismatch')
        candidates = [r for r in pr['batches'] if r['branch'] == branch and r['stage'] == 'relax']
        require(len(candidates) == 1, 'ambiguous parent')
        row = candidates[0]
        require(row['cumulative_sweeps'] == 8 and row['unpinned_sweeps'] == (8 if branch == 'random' else 4)
                and row['settings']['maxdim'] == [128, 128], 'wrong parent stage')
        require(rootpath(row['checkpoint']) == rootpath(spec['checkpoint'])
                and row['checkpoint_sha256']['metadata.toml'] == spec['metadata_sha256'], 'wrong selected parent')
        meta, checks[f'{branch}_parent8'] = audit_checkpoint(row, pr, pp.parent, pins)
        saved = record['parents'][branch]
        require(rootpath(saved['record']) == pp and saved['record_sha256'] == spec['record_sha256']
                and rootpath(saved['checkpoint']) == rootpath(spec['checkpoint'])
                and saved['checkpoint_sha256'] == row['checkpoint_sha256']
                and saved['initial_completed_sweeps'] == 8 and saved['settings'] == row['settings']
                and abs(saved['loaded_energy'] - row['energy']) <= 1e-10
                and saved['strict_source_runtime_configuration_load'] is True, 'parent load provenance mismatch')
        parents[branch], metas[branch] = row, meta
        templates = pr['templates'] if templates is None else templates
        require(templates == pr['templates'], 'parent template mismatch')
    require(metas['random']['state']['site_indices'] == metas['windmill']['state']['site_indices'],
            'parent physical site indices differ')
    children = {}
    for row in record['batches']:
        branch, chi = row['branch'], row['chi']
        require(row['cumulative_sweeps'] == 10 and row['unpinned_sweeps'] == (10 if branch == 'random' else 6),
                'wrong child sweep count')
        require(rootpath(row['parent_checkpoint']) == rootpath(parents[branch]['checkpoint'])
                and row['parent_record_sha256'] == record['parents'][branch]['record_sha256'],
                'child not tied to selected eight-sweep parent')
        expected_settings = {**parents[branch]['settings'], 'maxdim': [chi, chi]}
        require(row['settings'] == expected_settings, 'solver changes beyond declared maxdim')
        require(row['same_chi_comparison'] == (chi == 128)
                and 0 <= row['parent_profile_unchanged_error'] <= 1e-12,
                'matched-parent solver provenance failed')
        meta, checks[f'{branch}_chi{chi}_child10'] = audit_checkpoint(row, record, path.parent, pins)
        require(meta['state']['site_indices'] == metas[branch]['state']['site_indices'], 'child site indices changed')
        precision(parents[branch], row)
        children[(branch, chi)] = row
    return record, parents, children, templates, pins, checks, execution


def pair_comparison(geometry, left, right, edge):
    result = base.pair_comparison(geometry, left, right, edge)
    sites, _ = base.window(geometry, edge)
    pairs = [(i, j) for i in sites for j in sites if i != j]
    for entry in result['all_circumference_shifts']:
        mapping, _ = base.shift_permutations(geometry, entry['y_shift'])
        a, b = left['correlations'], right['correlations']
        components = {key: base.stats([b[key][mapping[i]][mapping[j]] - a[key][i][j]
                                      for i, j in pairs]) for key in ('zz_real', 'zz_imag', 'pm_real', 'pm_imag')}
        components['connected_zz_real'] = base.stats([
            b['zz_real'][mapping[i]][mapping[j]] - right['sz_profile'][mapping[i]] * right['sz_profile'][mapping[j]]
            - a['zz_real'][i][j] + left['sz_profile'][i] * left['sz_profile'][j] for i, j in pairs])
        for prefix in ('zz', 'pm'):
            components[prefix + '_complex_rms'] = math.hypot(components[prefix + '_real']['rms'],
                                                            components[prefix + '_imag']['rms'])
        entry['correlations'] = components
    result.update(edge_columns_excluded=edge, included_columns=list(range(edge, geometry['Lx']-edge)),
        correlation_selection='ordered off-site pairs; both endpoints in retained columns',
        correlation_alignment='both matrix indices translated by the same bond-RMS-selected shift',
        correlation_pair_count=len(pairs))
    return result


def compare(geometry, left, right):
    return dict(energy_change=right['energy'] - left['energy'],
                windows={str(e): pair_comparison(geometry, left, right, e) for e in EDGES})


def row_summary(geometry, row, templates):
    keys = ('energy', 'energy_per_site', 'variance', 'variance_per_site', 'maxlinkdim',
            'sz_profile', 'bond_energy', 'cumulative_sweeps', 'unpinned_sweeps', 'checkpoint')
    return dict(**{key: row[key] for key in keys},
        measured_truncation_error=row['measured_truncation_errors'][-1],
        windows={str(e): base.profile_summary(geometry, row, templates, e) for e in EDGES},
        correlation_lines=base.correlation_lines(geometry, row))


def analyze(record, parents, children, templates):
    g = record['configuration']
    result = dict(schema_version=1, scope='N54_Q6_equal_two_sweep_matched_parent_structure_comparison',
        model='nearest-neighbor isotropic J=1; theta=0; preparation lambda=0', geometry=g,
        selection='random/windmill selected by aligned central bond contrast, not trial energy',
        variance_interpretation='Hamiltonian dispersion, not a ground-energy error bar',
        chi_interpretation='same-parent finite two-sweep response; not a converged chi dependence',
        phase_identification='not_attempted', axial_alignment='not_tested', branches={}, branch_pairs={})
    for branch in BRANCHES:
        result['branches'][branch] = dict(parent8=row_summary(g, parents[branch], templates), children={},
            parent8_to_child10={}, chi128_to_chi256=compare(g, children[(branch, 128)], children[(branch, 256)]))
        for chi in CHIS:
            child = children[(branch, chi)]
            result['branches'][branch]['children'][str(chi)] = row_summary(g, child, templates)
            result['branches'][branch]['parent8_to_child10'][str(chi)] = dict(
                **compare(g, parents[branch], child), precision=precision(parents[branch], child),
                all_precision_conditions_passed=child['all_precision_conditions_passed'])
    result['branch_pairs']['parent8'] = compare(g, parents['random'], parents['windmill'])
    for chi in CHIS:
        result['branch_pairs'][str(chi)] = compare(g, children[('random', chi)], children[('windmill', chi)])
    # Compare separations with each branch's motion, using the same observable
    # support. These are descriptive ratios; motion is not an uncertainty bar.
    result['separation_relative_to_parent_child_motion'] = {}
    for chi in CHIS:
        entries = {}
        for edge in EDGES:
            pair = result['branch_pairs'][str(chi)]['windows'][str(edge)]
            comparisons = {}
            for alignment in ('raw', 'minimum_bond_rms_shift'):
                gap = pair[alignment]
                motion = {b: result['branches'][b]['parent8_to_child10'][str(chi)]['windows'][str(edge)][alignment]
                          for b in BRANCHES}
                comparisons[alignment] = {}
                for observable in ('bond', 'sz'):
                    diffs = {b: motion[b][observable]['rms'] for b in BRANCHES}
                    scale = max(diffs.values())
                    comparisons[alignment][observable] = dict(branch_rms=gap[observable]['rms'],
                        parent_child_rms=diffs, ratio_to_largest_motion=gap[observable]['rms']/scale if scale else None)
                for observable in ('connected_zz_real', 'pm_real', 'pm_imag'):
                    diffs = {b: motion[b]['correlations'][observable]['rms'] for b in BRANCHES}
                    scale = max(diffs.values())
                    comparisons[alignment][observable] = dict(branch_rms=gap['correlations'][observable]['rms'],
                        parent_child_rms=diffs, ratio_to_largest_motion=gap['correlations'][observable]['rms']/scale if scale else None)
            entries[str(edge)] = comparisons
        result['separation_relative_to_parent_child_motion'][str(chi)] = entries
    result['motion_ratio_interpretation'] = ('Descriptive scales, not statistical errors or phase criteria; '
        'each comparison selects its own bond-minimizing circumference shift, shared by all its observables.')
    passed = all(row['all_precision_conditions_passed'] for row in children.values())
    result['claim_status'] = ('selected_precision_checks_passed_chi_size_seed_convergence_unestablished' if passed
                              else 'unconverged_trials_no_energy_ranking_claim')
    return result


def render(result, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout='constrained')
    labels = ['parent 8', 'chi128 +2', 'chi256 +2']
    for col, edge in enumerate(EDGES):
        for axis, observable in zip(axes[:, col], ('bond', 'sz')):
            for alignment, style, name in [('raw', 'o--', 'branch difference, raw'),
                                          ('minimum_bond_rms_shift', 'o-', 'branch difference, aligned')]:
                axis.plot(range(3), [result['branch_pairs'][stage]['windows'][str(edge)][alignment][observable]['rms']
                                    for stage in ('parent8', '128', '256')], style, label=name)
            for branch, marker in zip(BRANCHES, ('s:', '^:')):
                axis.plot([1, 2], [result['branches'][branch]['parent8_to_child10'][str(chi)]['windows'][str(edge)]
                                  ['minimum_bond_rms_shift'][observable]['rms'] for chi in CHIS],
                          marker, label=f'{branch}: aligned parent-child motion')
            axis.set_xticks(range(3), labels)
            axis.set_ylabel(f'{observable} RMS difference')
            axis.set_title(f'Exclude {edge} edge column(s) at each end')
            axis.grid(alpha=.2)
            axis.legend(fontsize=8)
    fig.suptitle('N54 Q6 matched 8-sweep parents, equal additional 2 sweeps\n'
                 'Finite-sweep structure response; motion is not an error bar; no phase or energy ranking')
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    require(not path.exists(), f'refusing existing figure: {path}')
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return dict(matplotlib_version=matplotlib.__version__, figure_sha256={label(path): digest(path)})


def render_profiles(result, path):
    """Show all children on common scales, aligning each pair by its inner window."""
    import matplotlib.pyplot as plt
    g = result['geometry']
    descriptors = [base.bond_descriptor(g, bond) for bond in g['bonds']]
    kinds = sorted({d[0] for d in descriptors})
    rows = []
    for chi in CHIS:
        shift = result['branch_pairs'][str(chi)]['windows']['2']['minimum_bond_rms_shift']['y_shift']
        for branch in BRANCHES:
            row = result['branches'][branch]['children'][str(chi)]
            dy = shift if branch == 'windmill' else 0
            sm, bm = base.shift_permutations(g, dy)
            rows.append((branch, chi, dy, [row['sz_profile'][i] for i in sm],
                         [row['bond_energy'][i] for i in bm]))
    szlim = max(abs(v) for row in rows for v in row[3])
    blo, bhi = min(v for row in rows for v in row[4]), max(v for row in rows for v in row[4])
    fig, axes = plt.subplots(2, 4, figsize=(15, 8), layout='constrained')
    for col, (branch, chi, dy, sz, bonds) in enumerate(rows):
        szgrid = [[math.nan]*6 for _ in range(9)]
        for value, site in zip(sz, g['sites']):
            szgrid[3*site['y']+'ABC'.index(site['sublattice'])][site['x']] = value
        im = axes[0, col].imshow(szgrid, origin='lower', aspect='auto', cmap='RdBu_r', vmin=-szlim, vmax=szlim)
        axes[0, col].set_title(f'{branch}, chi={chi}, dy={dy}')
        axes[0, col].set_yticks(range(9), [f'{s}, y={y}' for y in range(3) for s in 'ABC'])
        grid = [[math.nan]*6 for _ in range(3*len(kinds))]
        for value, (kind, x, y) in zip(bonds, descriptors):
            grid[3*kinds.index(kind)+y][x] = value
        bm = axes[1, col].imshow(grid, origin='lower', aspect='auto', cmap='viridis', vmin=blo, vmax=bhi)
        axes[1, col].set_yticks([3*i+1 for i in range(len(kinds))], kinds, fontsize=8)
        for axis in axes[:, col]:
            axis.set_xticks(range(6))
            axis.axvline(.5, color='white', ls='--', lw=.9)
            axis.axvline(4.5, color='white', ls='--', lw=.9)
            axis.axvspan(1.5, 3.5, facecolor='none', edgecolor='#f4a340', lw=1.5)
        axes[0, col].set_xlabel('site column x')
        axes[1, col].set_xlabel('bond anchor column x')
    fig.colorbar(im, ax=axes[0, :], label='physical Sz (common scale)', shrink=.85)
    fig.colorbar(bm, ax=axes[1, :], label='bond energy / J (common scale)', shrink=.85)
    fig.suptitle('N54 Q6, nearest-neighbor trials after matched +2 sweeps\n'
        'Each chi pair aligned by bond RMS in x=2,3; raw data retained separately\n'
        'Orange: site/bond-anchor columns x=2,3; dashed: exclude one edge column\n'
        'Numerical bond windows require both endpoints inside the retained columns', fontsize=12)
    path = Path(path)
    require(not path.exists(), f'refusing existing figure: {path}')
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return {label(path): digest(path)}


def self_test():
    base.shift_self_test()
    sites = [dict(index=1+9*x+3*y+s, x=x, y=y, sublattice='ABC'[s])
             for x in range(6) for y in range(3) for s in range(3)]
    lookup = {(s['x'], s['y'], s['sublattice']): s['index'] for s in sites}
    bonds = []
    for x in range(6):
        for y in range(3):
            for a, b, dx, dy in [('A','B',0,0), ('A','C',0,0), ('B','C',0,0),
                                  ('A','B',-1,0), ('A','C',0,-1), ('B','C',1,-1)]:
                if 0 <= x+dx < 6:
                    bonds.append(dict(i=lookup[(x,y,a)], j=lookup[(x+dx,(y+dy)%3,b)], wy=(y+dy)//3))
    g = dict(Lx=6, Ly=3, N=54, sites=sites, bonds=bonds)
    left = dict(sz_profile=[i/100 for i in range(54)], bond_energy=[float(i*i) for i in range(102)],
        correlations={k: [[(i*i+7*j+3*i*j+c)/10000 for j in range(54)] for i in range(54)]
                      for c, k in enumerate(('zz_real','zz_imag','pm_real','pm_imag'))})
    sm, bm = base.shift_permutations(g, 1)
    right = dict(sz_profile=[0.0]*54, bond_energy=[0.0]*102,
                 correlations={k: [[0.0]*54 for _ in range(54)] for k in left['correlations']})
    for i, j in enumerate(sm):
        right['sz_profile'][j] = left['sz_profile'][i]
        for k, l in enumerate(sm):
            for name in left['correlations']:
                right['correlations'][name][j][l] = left['correlations'][name][i][k]
    for i, j in enumerate(bm):
        right['bond_energy'][j] = left['bond_energy'][i]
    for edge, count in ((1, 36*35), (2, 18*17)):
        result = pair_comparison(g, left, right, edge)
        best = result['minimum_bond_rms_shift']
        require(best['y_shift'] == 1 and best['bond']['rms'] == best['sz']['rms'] == 0,
                'synthetic profile alignment failed')
        require(result['correlation_pair_count'] == count
                and all(v['max_abs'] <= 1e-15 for v in best['correlations'].values() if isinstance(v, dict)),
                'two-index correlation translation failed')
        require(result['raw']['correlations']['pm_complex_rms'] > 0, 'fixture cannot expose missing translation')
    for matrix in right['correlations'].values():
        for i in range(54):
            matrix[i][i] += 100
    require(pair_comparison(g, left, right, 2)['minimum_bond_rms_shift']['correlations']['pm_complex_rms'] == 0,
            'on-site correlation terms were not excluded')
    return dict(status='passed', checks=['independent 102-bond winding fixture', 'both correlation indices translated',
        'common bond-selected shift for all observables', 'both edge windows', 'on-site correlation exclusion'],
        not_a_DMRG_result=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('validation', type=Path, nargs='?')
    parser.add_argument('--output', type=Path, default=ROOT/'docs/research/data/p4_vbc54_matched.json')
    parser.add_argument('--figure', type=Path, default=ROOT/'docs/research/figures/p4_vbc54_matched.png')
    parser.add_argument('--no-plot', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(), indent=2))
        if args.validation is None:
            return
    if args.validation is None:
        parser.error('validation is required unless --self-test')
    require(not args.output.exists(), f'refusing existing output: {args.output}')
    if not args.no_plot:
        require(not args.figure.exists(), f'refusing existing figure: {args.figure}')
        profiles = args.figure.with_name(args.figure.stem + '_profiles.png')
        require(not profiles.exists(), f'refusing existing figure: {profiles}')
    sources = {label(p): digest(p) for p in (Path(__file__), Path(base.__file__))}
    record, parents, children, templates, pins, checks, execution = load_audited(args.validation)
    result = analyze(record, parents, children, templates)
    result['provenance'] = dict(input_sha256=pins, analysis_source_sha256=sources, python_version=sys.version.split()[0],
        checks=checks, translation_fixture=self_test(), checkpoint_payloads_hash_verified=True, mps_recontracted=False,
        selected_parent_checkpoints_only=True, sources_unchanged=record['sources_unchanged'],
        backend_unchanged=record['backend_unchanged'], parents_unchanged=record['parents_unchanged'],
        execution={k: execution[k] for k in ('status', 'wall_elapsed_seconds', 'wall_limit_seconds',
                                            'worker_exit_code', 'worker_exit_confirmed', 'julia_threads', 'blas_threads')})
    if not args.no_plot:
        result['artifacts'] = render(result, args.figure)
        result['artifacts']['figure_sha256'].update(render_profiles(result, profiles))
    require(all(digest(rootpath(p)) == sha for p, sha in {**pins, **sources}.items()), 'input/source changed during analysis')
    result['provenance']['inputs_unchanged_after_analysis'] = True
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(result, out, ensure_ascii=False, indent=2, allow_nan=False)
        out.write('\n')
    print(json.dumps(dict(output=label(args.output), status=result['claim_status'],
        precision={b: {str(c): result['branches'][b]['parent8_to_child10'][str(c)]['precision']
                       for c in CHIS} for b in BRANCHES}), indent=2))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError) as error:
        raise SystemExit(str(error)) from error
