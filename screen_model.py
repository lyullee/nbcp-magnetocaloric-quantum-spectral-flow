"""Finite-cluster screening, not a bulk-material or quantum-advantage result.

Hamiltonian in kelvin (H/k_B):
sum_edges [Jxy*(Sx_i*Sx_j+Sy_i*Sy_j)+Jz*Sz_i*Sz_j]
- g*mu_B/k_B*B*sum_i Sz_i, with S=Pauli/2.

Input: Xiang et al., Nature 625,270 (2024), author manuscript Methods.
The full longitudinal-field grid reuses zero-field spectra because total Sz
commutes with the exchange Hamiltonian. Dense full-state diagonalization is
not needed: diagonalize magnetization sectors once.
"""
from collections import deque
from itertools import combinations
from pathlib import Path
import csv
import json
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
JXY = 0.88
JZ = 1.48
G = 4.89
MU_B_OVER_KB = 0.67171381563
R_GAS = 8.31446261815324
ION_DENSITY_PER_CM3 = 5.8e21
N_A = 6.02214076e23
VOLUME_FACTOR_MJ_PER_CM3 = 1e3 * ION_DENSITY_PER_CM3 / N_A * R_GAS

def lattice(a, b):
    """Triangular periodic quotient with exact integer keys, no wrap guessing."""
    determinant = a[0]*b[1]-a[1]*b[0]
    assert determinant > 0
    def key(x, y):
        return ((b[1]*x-b[0]*y) % determinant,
                (-a[1]*x+a[0]*y) % determinant)
    nodes = {key(0, 0): (0, 0)}
    queue = deque([(0, 0)])
    directions = [(1, 0), (0, 1), (1, -1)]
    while queue:
        x, y = queue.popleft()
        for dx, dy in directions:
            point = (x+dx, y+dy)
            k = key(*point)
            if k not in nodes:
                nodes[k] = point
                queue.append(point)
    assert len(nodes) == determinant
    indices = {k: i for i, k in enumerate(nodes)}
    edges = set()
    for k, (x, y) in nodes.items():
        for dx, dy in directions:
            pair = tuple(sorted((indices[k], indices[key(x+dx, y+dy)])))
            assert pair[0] != pair[1]
            edges.add(pair)
    edges = sorted(edges)
    degrees = [sum(i in edge for edge in edges) for i in range(determinant)]
    assert degrees == [6]*determinant and len(edges) == 3*determinant
    # Both geometries support the three-sublattice order: no artificial seam.
    colors = [(x+2*y) % 3 for x, y in nodes.values()]
    assert all(colors[i] != colors[j] for i, j in edges)
    return {"N": determinant, "translations": [a, b], "edges": edges,
            "coordinates": list(nodes.values()), "degree": degrees}

def spectra(graph):
    n = graph["N"]
    xxz_energy, mz, ising_energy = [], [], []
    max_dim = 0
    for k in range(n+1):
        basis = [sum(1 << i for i in occupied) for occupied in combinations(range(n), k)]
        index = {s: p for p, s in enumerate(basis)}
        max_dim = max(max_dim, len(basis))
        h = np.zeros((len(basis), len(basis)))
        for p, s in enumerate(basis):
            for i, j in graph["edges"]:
                equal = ((s >> i) & 1) == ((s >> j) & 1)
                h[p, p] += JZ * (0.25 if equal else -0.25)
                if not equal:
                    h[p, index[s ^ (1 << i) ^ (1 << j)]] += JXY / 2
        assert np.max(np.abs(h-h.T)) < 1e-12
        xxz_energy.extend(np.linalg.eigvalsh(h).tolist())
        ising_energy.extend(np.diag(h).tolist())
        mz.extend([k-n/2]*len(basis))
    graph["largest_sector_dimension"] = max_dim
    assert len(xxz_energy) == 2**n
    result = {"XXZ": np.asarray(xxz_energy), "Ising": np.asarray(ising_energy),
              "Mz": np.asarray(mz)}
    np.savez_compressed(OUT / f"spectrum_N{n}.npz", **result)
    return result

def thermal(energy_0, magnetization, n, temperature, fields):
    fields = np.atleast_1d(fields)
    energies = energy_0[:, None] - (G*MU_B_OVER_KB)*magnetization[:, None]*fields[None, :]
    minimum = energies.min(axis=0)
    w = np.exp(-(energies-minimum)/temperature)
    z = w.sum(axis=0)
    probability = w/z
    u = (probability*energies).sum(axis=0)
    entropy = (np.log(z)+(u-minimum)/temperature)/n
    cv = ((probability*(energies-u)**2).sum(axis=0)/temperature**2)/n
    m = (probability*magnetization[:, None]).sum(axis=0)/n
    return {"S_per_spin_kB": entropy, "C_per_spin_kB": cv,
            "U_per_spin_K": u/n, "Mz_per_spin": m}

def numeric_pair(frame, x, y):
    d = frame.iloc[5:, [x, y]].apply(pd.to_numeric, errors="coerce").dropna()
    return d.to_numpy(dtype=float)

def observed_audit():
    audit = {"source_workbooks": {}, "cooling_curves": [], "entropy_curves": [],
             "notes": ["Corrected points use experimental entropy and are not independent thermometer observations.",
                       "Theory columns are published XTRG calculations, not observations.",
                       "All spreadsheet headers/units have been inspected; original workbook bytes are unchanged."]}
    for f in (ROOT / "sources").glob("*.xlsx"):
        audit["source_workbooks"][f.name] = {
            name: {"rows": len(d), "columns": len(d.columns)}
            for name, d in pd.read_excel(f, sheet_name=None, header=None).items()
        }
    f = pd.read_excel(ROOT / "sources/figure_2.xlsx", sheet_name="Figure 2a", header=None)
    for x, y in [(0,1),(2,3),(4,5),(7,8),(11,12)]:
        a = numeric_pair(f, x, y)
        audit["cooling_curves"].append({"label": str(f.iloc[4,y]), "points": len(a),
            "minimum_T_K": float(a[:,1].min()), "field_at_minimum_T": float(a[a[:,1].argmin(),0])})
    f = pd.read_excel(ROOT / "sources/extended_figure_6.xlsx", sheet_name="c", header=None)
    entropy0 = numeric_pair(f,0,1)
    entropy4 = numeric_pair(f,2,3)
    s_hot = float(np.interp(2., entropy4[np.argsort(entropy4[:,0]),0], entropy4[np.argsort(entropy4[:,0]),1]))
    audit["measured_entropy_at_2K_4T_mJ_per_cm3_K"] = s_hot
    for x,y in [(0,1),(2,3)]:
        a = numeric_pair(f,x,y)
        audit["entropy_curves"].append({"label": str(f.iloc[4,y]), "points": len(a),
            "T_range_K": [float(a[:,0].min()),float(a[:,0].max())],
            "S_range_mJ_per_cm3_K": [float(a[:,1].min()),float(a[:,1].max())]})
    # Endpoint entropy budget: Q = T*(S_final - S_hot), conditional on a
    # reversible route reaching the cold hold leg. Not a full cycle simulation.
    for temp in [0.1, 0.2, 0.3]:
        a = entropy0[np.argsort(entropy0[:,0])]
        assert a[0,0] <= temp <= a[-1,0]
        s_end = float(np.interp(temp, a[:,0], a[:,1]))
        audit.setdefault("source_equilibrium_hold_leg", []).append({
            "T_cold_K": temp, "B_end_T": 0., "T_hot_K": 2., "B_hot_T": 4.,
            "Q_cold_mJ_per_cm3": temp*(s_end-s_hot), "S_end_mJ_per_cm3_K": s_end,
            "warning": "Endpoint entropy budget only. Reachability of the full cycle, heat switches, ramp time, phonons and nuclear heat are not modelled."})
    candidates = [{"B_end_T":0., "S_cold_mJ_per_cm3_K":float(np.interp(.1, entropy0[np.argsort(entropy0[:,0]),0],entropy0[np.argsort(entropy0[:,0]),1]))}]
    source = pd.read_excel(ROOT / "sources/extended_figure_4.xlsx",sheet_name=None,header=None)
    for sheet,b in [("h",1.3),("i",1.5),("j",1.66)]:
        a = numeric_pair(source[sheet],0,3)
        a = a[np.argsort(a[:,0])]
        assert a[0,0] <= .1 <= a[-1,0]
        candidates.append({"B_end_T":b, "S_cold_mJ_per_cm3_K":float(np.interp(.1,a[:,0],a[:,1]))})
    for c in candidates:
        c["T_cold_K"] = .1
        c["Q_entropy_budget_mJ_per_cm3"] = .1*(c["S_cold_mJ_per_cm3_K"]-s_hot)
    ranked = sorted(candidates,key=lambda c:c["Q_entropy_budget_mJ_per_cm3"],reverse=True)
    gap = ranked[0]["Q_entropy_budget_mJ_per_cm3"]-ranked[1]["Q_entropy_budget_mJ_per_cm3"]
    audit["measured_field_candidates_100mK"] = candidates
    audit["conditional_decision_accuracy_target"] = {
        "best_of_measured_candidates_T":ranked[0]["B_end_T"],
        "runner_up_of_measured_candidates_T":ranked[1]["B_end_T"],
        "Q_gap_mJ_per_cm3":gap,
        "equal_absolute_entropy_error_bound_each_endpoint_mJ_per_cm3_K":gap/(2*.1),
        "equal_absolute_entropy_error_bound_each_endpoint_kB_per_spin":gap/(2*.1)/VOLUME_FACTOR_MJ_PER_CM3,
        "note":"Sufficient central-value ranking bound for two endpoint entropies with shared hot state; not an experimental confidence interval and not a full-cycle optimum. Other input and interpolation uncertainties must be added."
    }
    return audit

def evaluate(n, data):
    fields = np.linspace(0., 2., 401)
    rows, profiles, checks, heat_capacity_comparison, entropy_comparison = [], [], {}, [], []
    for model in ["XXZ", "Ising"]:
        energy = data[model]
        m = data["Mz"]
        # Independent physical checks rather than self-comparison.
        q = thermal(energy, m, n, 1e6, [0.])
        checks[model] = {"high_T_entropy_error": float(abs(q["S_per_spin_kB"][0]-np.log(2))),
                         "high_T_energy_K_per_spin": float(q["U_per_spin_K"][0]),
                         "polarized_4T_Mz_per_spin_at_1mK": float(thermal(energy,m,n,.001,[4.])["Mz_per_spin"][0])}
        assert checks[model]["high_T_entropy_error"] < 1e-10
        assert abs(checks[model]["polarized_4T_Mz_per_spin_at_1mK"]-.5) < 1e-10
        t,b = .25,.31
        dt = 1e-5
        here = thermal(energy,m,n,t,[b])
        above = thermal(energy,m,n,t+dt,[b])
        below = thermal(energy,m,n,t-dt,[b])
        c_error = abs((above["U_per_spin_K"][0]-below["U_per_spin_K"][0])/(2*dt)-here["C_per_spin_kB"][0])
        s_error = abs((above["S_per_spin_kB"][0]-below["S_per_spin_kB"][0])/(2*dt)-here["C_per_spin_kB"][0]/t)
        checks[model].update({"dU_dT_minus_C":float(c_error), "dS_dT_minus_C_over_T":float(s_error)})
        assert c_error < 1e-6 and s_error < 1e-6
        s_hot = thermal(energy,m,n,2.,[4.])["S_per_spin_kB"][0]
        for temp in [.1,.2,.3,.5]:
            quantities = thermal(energy,m,n,temp,fields)
            s = quantities["S_per_spin_kB"]
            q_cold = temp*(s-s_hot)*VOLUME_FACTOR_MJ_PER_CM3
            i = int(np.argmax(q_cold))
            end_energies = energy-(G*MU_B_OVER_KB)*m*fields[i]
            residual = np.log(np.sum(np.abs(end_energies-end_energies.min())<1e-9))/n
            rows.append({"N":n,"model":model,"T_cold_K":temp,"T_hot_K":2.,"B_hot_T":4.,
                         "B_end_opt_grid_T":float(fields[i]), "Q_opt_mJ_per_cm3":float(q_cold[i]),
                         "Q_B0_mJ_per_cm3":float(q_cold[0]), "S_hot_kB_per_spin":float(s_hot),
                         "end_ground_state_entropy_kB_per_spin":float(residual),
                         "ground_entropy_exceeds_hot_entropy":bool(residual>s_hot),
                         "scope":"Endpoint entropy budget, not verified cycle reachability"})
            for j,b_field in enumerate(fields):
                profiles.append({"N":n,"model":model,"T_cold_K":temp,"B_end_T":float(b_field),
                                 "Q_entropy_budget_mJ_per_cm3":float(q_cold[j]),
                                 "S_kB_per_spin":float(s[j])})
        # Compare small-cluster thermodynamics with published experimental heat capacity.
        source = pd.read_excel(ROOT / "sources/figure_2.xlsx", sheet_name="Figure 2c", header=None)
        for x,y,b_field in [(0,1,0.),(3,4,.33),(6,7,.75),(9,10,1.66)]:
            for t_measured, cm_over_t in numeric_pair(source,x,y):
                if .06 <= t_measured <= 2.:
                    predicted = thermal(energy,m,n,float(t_measured),[b_field])["C_per_spin_kB"][0]*R_GAS/t_measured
                    heat_capacity_comparison.append({"N":n,"model":model,"B_T":b_field,"T_K":float(t_measured),
                        "source_Cm_over_T_J_per_mol_K2":float(cm_over_t),"finite_cluster_Cm_over_T_J_per_mol_K2":float(predicted)})
        source = pd.read_excel(ROOT / "sources/extended_figure_4.xlsx", sheet_name=None, header=None)
        for sheet,b_field in [("h",1.3),("i",1.5),("j",1.66)]:
            for t_measured, s_measured in numeric_pair(source[sheet],0,3):
                if .06 <= t_measured <= 2.:
                    s_predicted = thermal(energy,m,n,float(t_measured),[b_field])["S_per_spin_kB"][0]*VOLUME_FACTOR_MJ_PER_CM3
                    entropy_comparison.append({"N":n,"model":model,"B_T":b_field,"T_K":float(t_measured),
                        "source_S_mJ_per_cm3_K":float(s_measured),"finite_cluster_S_mJ_per_cm3_K":float(s_predicted)})
    return rows,profiles,checks,heat_capacity_comparison,entropy_comparison

def save_csv(name, rows):
    with (OUT / name).open("w",newline="",encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def main():
    start = time.perf_counter()
    graphs = [lattice((3,0),(0,3)), lattice((2,2),(-2,4))]
    all_rows, all_profiles, all_heat_capacity, all_entropy = [],[],[],[]
    checks = {}
    for graph in graphs:
        n = graph["N"]
        print(f"Computing N={n}, {len(graph['edges'])} bonds",flush=True)
        d = spectra(graph)
        rows,profiles,check,heat_capacity,entropy = evaluate(n,d)
        checks[str(n)] = check
        all_rows.extend(rows); all_profiles.extend(profiles)
        all_heat_capacity.extend(heat_capacity); all_entropy.extend(entropy)
        print(f"N={n} complete, elapsed={time.perf_counter()-start:.1f}s",flush=True)
    save_csv("finite_cluster_cycle_summary.csv",all_rows)
    save_csv("finite_cluster_cycle_profiles.csv",all_profiles)
    save_csv("heat_capacity_comparison.csv",all_heat_capacity)
    save_csv("entropy_comparison.csv",all_entropy)
    audit = observed_audit()
    summary = {"date":"2026-10-11","parameters":{"Jxy_K":JXY,"Jz_K":JZ,"g_c":G,
        "mu_B_over_kB_K_per_T":MU_B_OVER_KB,"magnetic_ions_per_cm3":ION_DENSITY_PER_CM3},
        "graphs":graphs,"physical_checks":checks,"source_data_audit":audit,
        "cycle_screen":all_rows,"elapsed_seconds":time.perf_counter()-start,
        "scope":["Exact diagonalization on a classical PC, not a quantum circuit or hardware execution.",
                 "N=9 and N=12 finite clusters are exact for those graphs, not a bulk convergence demonstration.",
                 "Ising is an ablation with Jxy set to zero, not the strongest classical algorithm for XXZ.",
                 "Q is the endpoint entropy budget conditional on a reversible cold hold leg. Full-cycle reachability has not been demonstrated; not device cooling power or electrical COP.",
                 "No new proof of quantum advantage and no finite-rate spin-phonon model.",
                 "Longitudinal-field sweeps commute with the XXZ exchange Hamiltonian; a closed-system ramp alone cannot establish thermalization."]}
    (OUT / "screening_summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"cycle_screen":all_rows,"source_data_audit":audit,
                      "elapsed_seconds":summary["elapsed_seconds"]},ensure_ascii=False,indent=2))

if __name__ == "__main__":
    main()
