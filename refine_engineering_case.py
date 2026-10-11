"""Public-data engineering benchmark, conditional magnetic-entropy budgets.

No field interpolation or temperature extrapolation. Candidate-field ranking
uses piecewise-linear interpolation of published central values. Scenario
loads are specified assumptions, not measured cooling-power data.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from screen_model import ROOT, OUT, numeric_pair, save_csv, VOLUME_FACTOR_MJ_PER_CM3

def prepared(values):
    order=np.argsort(values[:,0])
    values=values[order]
    temperatures, first, counts=np.unique(values[:,0],return_index=True,return_counts=True)
    # Published Fig4h repeats 12 identical temperature/entropy pairs. Collapse
    # only duplicates with agreeing entropy values; reject conflicting pairs.
    for t,count in zip(temperatures,counts):
        if count>1:
            assert np.ptp(values[values[:,0]==t,1])<1e-9
    values=values[first]
    assert np.all(np.diff(values[:,0])>0)
    return values

def interp(values,t):
    assert values[0,0]<=t<=values[-1,0]
    return float(np.interp(t,values[:,0],values[:,1]))

def roots_for_difference(a,b,lo,hi):
    knots=np.unique(np.concatenate(([lo,hi],a[:,0],b[:,0])))
    knots=knots[(knots>=lo)&(knots<=hi)]
    da=np.interp(knots,a[:,0],a[:,1])-np.interp(knots,b[:,0],b[:,1])
    roots=[]
    for i in range(len(knots)-1):
        if da[i]*da[i+1]<0:
            roots.append(float(knots[i]-da[i]*(knots[i+1]-knots[i])/(da[i+1]-da[i])))
    return roots

def roots_for_entropy(values,s):
    roots=[]
    for i in range(len(values)-1):
        x0,y0=values[i];x1,y1=values[i+1]
        if (y0-s)*(y1-s)<0:
            roots.append(float(x0+(s-y0)*(x1-x0)/(y1-y0)))
    return roots

def main():
    ef4=pd.read_excel(ROOT/'sources/extended_figure_4.xlsx',sheet_name=None,header=None)
    ef6=pd.read_excel(ROOT/'sources/extended_figure_6.xlsx',sheet_name='c',header=None)
    curves={0.:prepared(numeric_pair(ef6,0,1))}
    hot=prepared(numeric_pair(ef6,2,3))
    for sheet,b in [('h',1.3),('i',1.5),('j',1.66)]:
        curves[b]=prepared(numeric_pair(ef4[sheet],0,3))
    s_hot=interp(hot,2.)
    # Published bulk XTRG curves are a serious classical reference, kept
    # separate from experiments. Only three matched fields are available.
    theory={b:prepared(numeric_pair(ef4[sheet],7,9))
            for sheet,b in [('h',1.3),('i',1.5),('j',1.66)]}
    theory_hot=prepared(numeric_pair(ef4['h'],7,8))
    s_hot_theory=interp(theory_hot,2.)
    rows=[]
    ranking=[]
    for t in [.09,.1,.12,.15,.2,.25,.3]:
        at_t=[]
        for b,curve in curves.items():
            s=interp(curve,t)
            row={'T_cold_K':t,'B_end_T':b,'S_cold_mJ_per_cm3_K':s,
                 'Q_magnetic_entropy_budget_mJ_per_cm3':t*(s-s_hot)}
            rows.append(row);at_t.append(row)
        ordered=sorted(at_t,key=lambda r:r['Q_magnetic_entropy_budget_mJ_per_cm3'],reverse=True)
        gap=ordered[0]['Q_magnetic_entropy_budget_mJ_per_cm3']-ordered[1]['Q_magnetic_entropy_budget_mJ_per_cm3']
        ranking.append({'T_cold_K':t,'best_candidate_B_T':ordered[0]['B_end_T'],
            'runner_up_B_T':ordered[1]['B_end_T'],'gap_Q_mJ_per_cm3':gap,
            'conditional_equal_entropy_error_bound_kB_per_spin':gap/(2*t)/VOLUME_FACTOR_MJ_PER_CM3})
        qbest=ordered[0]['Q_magnetic_entropy_budget_mJ_per_cm3']
        ranking[-1]['candidates_within_5_percent_of_best_Q_T']=[r['B_end_T'] for r in ordered
            if r['Q_magnetic_entropy_budget_mJ_per_cm3']>=.95*qbest]
        ranking[-1]['runner_up_relative_Q_loss_percent']=100*gap/qbest
    save_csv('temperature_candidate_comparison.csv',rows)
    save_csv('temperature_candidate_ranking.csv',ranking)
    classical_reference=[]
    for t in [.09,.1,.12,.15,.2,.25,.3]:
        picked=max(theory,key=lambda b:interp(theory[b],t))
        observed_best=max(theory,key=lambda b:interp(curves[b],t))
        observed_q=t*(interp(curves[picked],t)-s_hot)
        observed_best_q=t*(interp(curves[observed_best],t)-s_hot)
        classical_reference.append({'T_cold_K':t,'published_XTRG_selected_B_T':picked,
            'experiment_based_selected_B_among_same_three_fields_T':observed_best,
            'observed_Q_budget_of_XTRG_choice_mJ_per_cm3':observed_q,
            'observed_best_Q_budget_mJ_per_cm3':observed_best_q,
            'choice_regret_Q_percent':100*(observed_best_q-observed_q)/observed_best_q,
            'scope':'Existing literature model vs experimental central values; NOT algorithm error of XTRG or quantum advantage.'})
    save_csv('published_XTRG_decision_comparison.csv',classical_reference)
    crossings=[]
    bs=list(curves)
    for i,b1 in enumerate(bs):
        for b2 in bs[i+1:]:
            for t in roots_for_difference(curves[b1],curves[b2],.09,.3):
                best=max(curves,key=lambda b:interp(curves[b],t))
                smax=interp(curves[best],t)
                crossings.append({'fields_T':[b1,b2],'T_K':t,
                    'both_on_central_value_best_candidate_envelope':abs(interp(curves[b1],t)-smax)<1e-8})
    scenarios=[]
    for row in rows:
        if row['T_cold_K']!=.1:continue
        q=row['Q_magnetic_entropy_budget_mJ_per_cm3']
        assert q>0
        heat_load_uW=.1
        hold_seconds=3600
        required_mJ=heat_load_uW*1e-3*hold_seconds
        scenarios.append({'T_cold_K':.1,'B_end_T':row['B_end_T'],
          'assumed_total_cold_load_uW':heat_load_uW,'required_hold_minutes':60,
          'required_heat_absorption_mJ':required_mJ,
          'conditional_magnetic_only_material_volume_cm3':required_mJ/q,
          'conditional_hold_minutes_for_0p4_cm3':q*.4/(heat_load_uW*1e-3)/60})
    save_csv('conditional_design_scenarios.csv',scenarios)
    reachability=[{'B_T':b,'T_bounds_K':[float(c[0,0]),float(c[-1,0])],
        'S_bounds_mJ_per_cm3_K':[float(c[:,1].min()),float(c[:,1].max())],
        'roots_S_equals_hot_T_K':roots_for_entropy(c,s_hot),
        'scope':'Magnetic-only isentropic endpoints within measured temperature range. Not full-cycle path validation.'}
        for b,c in curves.items()]
    result={'hot_state':{'T_K':2.,'B_T':4.,'S_mJ_per_cm3_K':s_hot},
        'comparisons':rows,'rankings':ranking,'candidate_crossings':crossings,
        'published_classical_reference_decisions':classical_reference,
        'design_scenarios':scenarios,'magnetic_isentrope_endpoint_checks':reachability,
        'scope':['Published central-value curves; no experimental confidence intervals assigned.',
                 'Twelve repeated temperature/entropy pairs in Fig4h agree within 1e-9; retained one copy in the derived curve. Source workbook unchanged.',
                 'No interpolation across magnetic fields. Candidate rankings only.',
                 'Magnetic entropy budgets only; nuclear, lattice, heat-switch and thermal-contact effects need separate modelling.',
                 'Device size and hold time are conditional design scenarios, not measured performance.',
                 'No quantum hardware or quantum circuit simulation.']}
    (OUT/'engineering_case_refinement.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ['comparisons']},indent=2))

if __name__=='__main__':main()
