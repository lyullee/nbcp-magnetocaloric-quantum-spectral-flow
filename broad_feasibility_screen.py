"""Broad finite-cluster feasibility screen for the proposed correction study.

This is a classical exact-diagonalization calculation.  It is deliberately a
screening calculation, not a prediction of the bulk material and not a quantum
hardware result.  The 2026 representative in-plane Hamiltonian is decomposed
as H = H_XXZ + a*V_Jpp + b*V_NN-dipole.  We scan a,b in [0,2].  Values other
than one are sensitivity factors, not fitted confidence intervals.

The engineering observable is the correction to a conditional magnetic
entropy heat budget between (2 K, 4 T) and (T_cold, B_end):

  delta Q = T_cold * [(S_full,cold-S_XXZ,cold)
                      -(S_full,hot-S_XXZ,hot)].

The public experimental heat budget is used only as a scale for judging
whether a correction is potentially material.  The two models use different
published parameter sets and are never added to make a hybrid prediction.
"""
import json
import time

import numpy as np

from screen_model import OUT, VOLUME_FACTOR_MJ_PER_CM3, lattice, save_csv
from screen_weak_terms import matrices, spectrum, entropy


TEMPERATURES = [.09, .10, .12, .15, .20, .25, .30]
FIELDS = [0., 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.66, 1.7, 1.8]
SCALES = [0., .5, 1., 1.5, 2.]
CHECK_FIELDS = [1.3, 1.5, 1.66]


def components(graph):
    """Return H0, V_Jpp, V_NN-dipole in the same parity bases."""
    base = matrices(graph, 'XXZ_2026_parameters')
    jpp = matrices(graph, 'XXZ_plus_additional_Jpp')
    full = matrices(graph, 'XXZ_plus_nn_dipoles_plus_additional_Jpp')
    out = []
    for (h0, mz0), (hj, mzj), (hf, mzf) in zip(base, jpp, full):
        assert np.array_equal(mz0, mzj) and np.array_equal(mz0, mzf)
        vj = hj-h0
        vd = hf-hj
        assert np.max(np.abs(vj-vj.conj().T)) < 1e-12
        assert np.max(np.abs(vd-vd.conj().T)) < 1e-12
        out.append((h0, vj, vd, mz0))
    return out


def assemble(parts, jpp_scale, dipole_scale):
    return [(h0+jpp_scale*vj+dipole_scale*vd, mz)
            for h0, vj, vd, mz in parts]


def public_scale_by_temperature():
    source = json.loads((OUT/'engineering_case_refinement.json').read_text(encoding='utf-8'))
    result = {}
    for t in TEMPERATURES:
        rows = [r for r in source['comparisons'] if abs(r['T_cold_K']-t) < 1e-12]
        assert rows
        result[t] = max(r['Q_magnetic_entropy_budget_mJ_per_cm3'] for r in rows)
        assert result[t] > 0
    return result


def evaluate_model(blocks, n, fields):
    hot_e = spectrum(blocks, 4.)
    hot_s = entropy(hot_e, 2., n)
    cold = {}
    for b in fields:
        e = spectrum(blocks, b)
        cold[b] = {t: entropy(e, t, n) for t in TEMPERATURES}
    return hot_s, cold


def main():
    start = time.perf_counter()
    public_scale = public_scale_by_temperature()
    graph9 = lattice((3,0),(0,3))
    parts9 = components(graph9)
    base9 = assemble(parts9, 0., 0.)
    base_hot9, base_cold9 = evaluate_model(base9, 9, FIELDS)

    rows = []
    for a in SCALES:
        for bscale in SCALES:
            blocks = assemble(parts9, a, bscale)
            hot, cold = evaluate_model(blocks, 9, FIELDS)
            for field in FIELDS:
                for t in TEMPERATURES:
                    delta_s_cold = cold[field][t]-base_cold9[field][t]
                    delta_s_hot = hot-base_hot9
                    delta_q = t*(delta_s_cold-delta_s_hot)*VOLUME_FACTOR_MJ_PER_CM3
                    base_q = t*(base_cold9[field][t]-base_hot9)*VOLUME_FACTOR_MJ_PER_CM3
                    design_scale = public_scale[t]
                    rows.append({
                        'N': 9,
                        'Jpp_scale': a,
                        'NN_dipole_scale': bscale,
                        'B_end_T': field,
                        'T_cold_K': t,
                        'delta_S_cold_kB_per_spin': delta_s_cold,
                        'delta_S_hot_kB_per_spin': delta_s_hot,
                        'base_XXZ_Q_model_mJ_per_cm3': base_q,
                        'perturbed_Q_model_mJ_per_cm3': base_q+delta_q,
                        'delta_Q_model_mJ_per_cm3': delta_q,
                        'public_best_Q_scale_mJ_per_cm3': design_scale,
                        'abs_delta_Q_as_percent_of_public_best_Q': 100*abs(delta_q)/design_scale,
                        'exceeds_10_percent_screen': bool(abs(delta_q) >= .1*design_scale),
                    })
            print(f'N=9 scales ({a:g},{bscale:g}) complete; elapsed={time.perf_counter()-start:.1f}s', flush=True)
    save_csv('broad_parameter_screen_N9.csv', rows)

    # Does re-optimizing the endpoint field remove the apparently large
    # fixed-field effect?  Compare each perturbed model with the base XXZ model
    # on exactly the same coarse field grid.
    decisions = []
    for a in SCALES:
        for bscale in SCALES:
            for t in TEMPERATURES:
                subset = [r for r in rows if r['Jpp_scale']==a and
                          r['NN_dipole_scale']==bscale and r['T_cold_K']==t]
                base_best = max(subset, key=lambda r:r['base_XXZ_Q_model_mJ_per_cm3'])
                full_best = max(subset, key=lambda r:r['perturbed_Q_model_mJ_per_cm3'])
                dq_opt = full_best['perturbed_Q_model_mJ_per_cm3']-base_best['base_XXZ_Q_model_mJ_per_cm3']
                decisions.append({
                    'N': 9,
                    'Jpp_scale': a,
                    'NN_dipole_scale': bscale,
                    'T_cold_K': t,
                    'base_optimal_B_on_grid_T': base_best['B_end_T'],
                    'perturbed_optimal_B_on_grid_T': full_best['B_end_T'],
                    'field_shift_T': full_best['B_end_T']-base_best['B_end_T'],
                    'base_optimal_Q_model_mJ_per_cm3': base_best['base_XXZ_Q_model_mJ_per_cm3'],
                    'perturbed_optimal_Q_model_mJ_per_cm3': full_best['perturbed_Q_model_mJ_per_cm3'],
                    'optimized_delta_Q_model_mJ_per_cm3': dq_opt,
                    'abs_optimized_delta_Q_as_percent_of_public_best_Q': 100*abs(dq_opt)/public_scale[t],
                    'optimized_effect_exceeds_10_percent_screen': bool(abs(dq_opt)>=.1*public_scale[t]),
                })
    save_csv('reoptimized_field_decisions_N9.csv', decisions)

    # Size check at the representative published scale.  This is intentionally
    # focused because N=12 dense diagonalization is the expensive part.
    size_rows = []
    for n, graph, parts in [(9, graph9, parts9),
                            (12, lattice((2,2),(-2,4)), None)]:
        if parts is None:
            parts = components(graph)
        base = assemble(parts, 0., 0.)
        realistic = assemble(parts, 1., 1.)
        base_hot, base_cold = evaluate_model(base, n, CHECK_FIELDS)
        full_hot, full_cold = evaluate_model(realistic, n, CHECK_FIELDS)
        for field in CHECK_FIELDS:
            for t in TEMPERATURES:
                delta_s_cold = full_cold[field][t]-base_cold[field][t]
                delta_s_hot = full_hot-base_hot
                delta_q = t*(delta_s_cold-delta_s_hot)*VOLUME_FACTOR_MJ_PER_CM3
                base_q = t*(base_cold[field][t]-base_hot)*VOLUME_FACTOR_MJ_PER_CM3
                size_rows.append({
                    'N': n,
                    'Jpp_scale': 1.,
                    'NN_dipole_scale': 1.,
                    'B_end_T': field,
                    'T_cold_K': t,
                    'delta_S_cold_kB_per_spin': delta_s_cold,
                    'delta_S_hot_kB_per_spin': delta_s_hot,
                    'base_XXZ_Q_model_mJ_per_cm3': base_q,
                    'perturbed_Q_model_mJ_per_cm3': base_q+delta_q,
                    'delta_Q_model_mJ_per_cm3': delta_q,
                    'abs_delta_Q_as_percent_of_public_best_Q': 100*abs(delta_q)/public_scale[t],
                })
        print(f'N={n} size check complete; elapsed={time.perf_counter()-start:.1f}s', flush=True)
    save_csv('published_scale_size_check.csv', size_rows)

    # Pair N=9 and N=12 to expose whether sign and rough magnitude survive.
    paired = []
    for field in CHECK_FIELDS:
        for t in TEMPERATURES:
            r9 = next(r for r in size_rows if r['N']==9 and r['B_end_T']==field and r['T_cold_K']==t)
            r12 = next(r for r in size_rows if r['N']==12 and r['B_end_T']==field and r['T_cold_K']==t)
            q9, q12 = r9['delta_Q_model_mJ_per_cm3'], r12['delta_Q_model_mJ_per_cm3']
            paired.append({
                'B_end_T': field,
                'T_cold_K': t,
                'delta_Q_N9_mJ_per_cm3': q9,
                'delta_Q_N12_mJ_per_cm3': q12,
                'same_sign': bool(np.sign(q9)==np.sign(q12)),
                'absolute_difference_mJ_per_cm3': abs(q12-q9),
                'N12_abs_percent_of_public_best_Q': 100*abs(q12)/public_scale[t],
                'N12_exceeds_10_percent_screen': bool(abs(q12)>=.1*public_scale[t]),
            })
    save_csv('finite_size_stability_check.csv', paired)

    realistic9 = [r for r in rows if r['Jpp_scale']==1. and r['NN_dipole_scale']==1.]
    nonzero9 = [r for r in rows if r['Jpp_scale']!=0. or r['NN_dipole_scale']!=0.]
    max_realistic9 = max(realistic9, key=lambda r:r['abs_delta_Q_as_percent_of_public_best_Q'])
    max_all9 = max(nonzero9, key=lambda r:r['abs_delta_Q_as_percent_of_public_best_Q'])
    realistic_decisions = [r for r in decisions if r['Jpp_scale']==1. and r['NN_dipole_scale']==1.]
    max_reoptimized_realistic = max(realistic_decisions,
        key=lambda r:r['abs_optimized_delta_Q_as_percent_of_public_best_Q'])
    n12_hits = [r for r in paired if r['N12_exceeds_10_percent_screen']]
    same_sign = sum(r['same_sign'] for r in paired)
    summary = {
        'calculation_type': 'Classical exact diagonalization on finite periodic clusters',
        'parameter_grid': {
            'temperatures_K': TEMPERATURES,
            'fields_T': FIELDS,
            'Jpp_scale': SCALES,
            'nearest_neighbor_dipole_scale': SCALES,
            'meaning_of_scale_one': 'Representative 2026 preprint parameter value',
            'meaning_of_other_scales': 'Sensitivity factors only; not statistical confidence bounds',
        },
        'engineering_metric': 'Correction to conditional magnetic entropy heat budget from (2 K,4 T) to (T,B)',
        'screening_threshold': 'abs(delta Q) >= 10% of the best public experimental central-value endpoint budget at the same T',
        'N9_grid_rows': len(rows),
        'N9_fraction_over_threshold': sum(r['exceeds_10_percent_screen'] for r in nonzero9)/len(nonzero9),
        'N9_representative_scale_maximum': max_realistic9,
        'N9_representative_scale_reoptimized_maximum': max_reoptimized_realistic,
        'N9_representative_scale_field_choices': realistic_decisions,
        'N9_entire_sensitivity_grid_maximum': max_all9,
        'N12_representative_points_over_threshold': len(n12_hits),
        'N12_representative_total_points': len(paired),
        'N9_N12_same_sign_points': same_sign,
        'N9_N12_total_points': len(paired),
        'N12_largest_representative_effect': max(paired, key=lambda r:r['N12_abs_percent_of_public_best_Q']),
        'scope': [
            'The public experimental budget supplies a comparison scale only; it is not added to the 2026 model correction.',
            'N=9 and N=12 results are exact for those clusters but are not thermodynamic-limit predictions.',
            'Only in-plane additional Jpp and nearest-neighbor dipoles are screened.',
            'Long-range dipoles, non-unique interlayer couplings, the incommensurate fan state, phonons, nuclear heat, heat switches and parasitic loads are absent.',
            'A threshold hit means the term deserves a stronger calculation; it does not validate a material or establish quantum advantage.',
            'No quantum circuit, quantum simulator or quantum hardware is used here.',
        ],
        'elapsed_seconds': time.perf_counter()-start,
    }
    (OUT/'broad_feasibility_screen.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main()
