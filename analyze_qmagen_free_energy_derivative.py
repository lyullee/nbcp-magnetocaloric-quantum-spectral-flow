"""Diagnose YC6-width entropy from staggered free-energy points near 90 mK."""
import json
from pathlib import Path

import numpy as np

from compile_qmagen_cylinder_convergence import public_reference

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def derivative_three_point(x: np.ndarray, y: np.ndarray, x0: float) -> float:
    """Derivative of the quadratic interpolant through three nonuniform points."""
    total = 0.0
    for i in range(3):
        others = [j for j in range(3) if j != i]
        numerator = (x0 - x[others[0]]) + (x0 - x[others[1]])
        denominator = (x[i] - x[others[0]]) * (x[i] - x[others[1]])
        total += y[i] * numerator / denominator
    return float(total)


def main(d: int = 24) -> None:
    paths = {
        0: RESULTS / f"qmagen_cylinder_9x6_D{d}.json",
        1: RESULTS / f"qmagen_cylinder_9x6_D{d}_st1.json",
        3: RESULTS / f"qmagen_cylinder_9x6_D{d}_st3.json",
    }
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing stagger runs: " + ", ".join(missing))
    data = {st: json.loads(path.read_text(encoding="utf-8")) for st, path in paths.items()}
    target = data[0]["cold_row"]
    beta0 = float(target["beta_dimensionless"])
    logz0 = float(target["logZ_per_spin"])

    all_points = []
    for st, run in data.items():
        for row in run["rows"]:
            all_points.append((float(row["beta_dimensionless"]), float(row["logZ_per_spin"]), st))
    all_points.sort()
    left = max(point for point in all_points if point[0] < beta0)
    right = min(point for point in all_points if point[0] > beta0)
    x = np.array([left[0], beta0, right[0]])
    y = np.array([left[1], logz0, right[1]])

    slope_left = (logz0 - left[1]) / (beta0 - left[0])
    slope_right = (right[1] - logz0) / (right[0] - beta0)
    slope_quad = derivative_three_point(x, y, beta0)
    entropy_left = logz0 - beta0 * slope_left
    entropy_right = logz0 - beta0 * slope_right
    entropy_quad = logz0 - beta0 * slope_quad
    convexity_consistent = bool(slope_left <= slope_right)
    direct = float(target["entropy_kB_per_spin"])
    ref = public_reference(float(target["temperature_K"]), 1.50)

    payload = {
        "method": "three staggered XTRG free-energy points around the st0 cold point",
        "bond_dimension": d,
        "target_temperature_K": float(target["temperature_K"]),
        "target_beta_dimensionless": beta0,
        "points": [
            {"beta": float(xx), "logZ_per_spin": float(yy), "stagger_index": int(st)}
            for xx, yy, st in (left, (beta0, logz0, 0), right)
        ],
        "slope_dlogZ_dbeta": {
            "left_secant": float(slope_left),
            "quadratic": float(slope_quad),
            "right_secant": float(slope_right),
        },
        "local_logZ_convexity_consistent": convexity_consistent,
        "entropy_kB_per_spin": {
            "direct_logZ_plus_betaE": direct,
            "free_energy_left_secant": float(entropy_left),
            "free_energy_quadratic": float(entropy_quad),
            "free_energy_right_secant": float(entropy_right),
            "secant_bracket_min": float(min(entropy_left, entropy_right)),
            "secant_bracket_max": float(max(entropy_left, entropy_right)),
            "public_reconstructed_proxy": ref,
        },
        "decision": (
            "DIAGNOSTIC_REJECTED_LOCAL_LOGZ_CONVEXITY_VIOLATION"
            if not convexity_consistent
            else "DIAGNOSTIC_ONLY_NOT_A_D_CONVERGENCE_TEST"
        ),
        "interpretation": (
            "For a physical partition function logZ(beta) is convex, so its secant slope must not decrease with beta. "
            "The independent stagger runs violate that condition; their truncation bias is larger than the local thermal curvature. "
            "The quadratic entropy estimate must therefore not be used as a result."
            if not convexity_consistent
            else "Diagnostic only: agreement of the free-energy derivative with the public proxy does not establish D convergence."
        ),
    }
    out = RESULTS / f"qmagen_cylinder_9x6_D{d}_free_energy_derivative.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
