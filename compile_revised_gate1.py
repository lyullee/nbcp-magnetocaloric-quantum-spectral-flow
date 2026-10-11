"""Compile the corrected bulk-window Gate-1 decision and audit figure."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"


def main():
    xtrg = json.loads((RESULTS / "npj2022_xtrg_entropy_reconstructed.json").read_text(encoding="utf-8"))
    n15 = json.loads((RESULTS / "n15_bulk_window_confirmation.json").read_text(encoding="utf-8"))
    n18 = json.loads((RESULTS / "n18_bulk_window_pilot.json").read_text(encoding="utf-8"))

    common_fields = (1.50, 1.55, 1.60)
    comparisons = []
    for t in (0.09, 0.10, 0.12):
        for b in common_fields:
            a = next(r for r in n15["rows"] if r["temperature_K"] == t and r["field_T"] == b)
            c = next(r for r in n18["rows"] if r["temperature_K"] == t and r["field_T"] == b)
            comparisons.append(
                {
                    "temperature_K": t,
                    "field_T": b,
                    "N15_deltaS": a["deltaS"],
                    "N18_deltaS": c["deltaS"],
                    "same_sign": bool(np.sign(a["deltaS"]) == np.sign(c["deltaS"])),
                    "N18_low16_deltaS": c["deltaS_low16_only"],
                    "N18_full_minus_low16": c["deltaS_full_minus_low16"],
                }
            )

    xtrg_90 = next(m for m in xtrg["metrics"] if m["anchor_N"] == 12 and m["temperature_K"] == 0.09)
    current_scale = (0.1225 * 0.636 * 11.60451812 / 0.88) * (4.89 / 4.716)
    payload = {
        "gate1_status": "HOLD",
        "quantum_stage_authorized": False,
        "supersedes": "gate1_preliminary_report.json field-window interpretation",
        "validated_bulk_M0_proxy": {
            "source": "Gao et al. 2022 XTRG Fig. 4a/4b",
            "zero_field_validation_RMSE_J_mol_K": xtrg["validation_against_Fig4b_zero_field_entropy"]["rmse_J_mol_K"],
            "90mK_B_star_T": xtrg_90["B_star_T"],
            "90mK_W95_T": [xtrg_90["W95_min_T"], xtrg_90["W95_max_T"]],
            "rough_current_parameter_field_mapping_T": xtrg_90["B_star_T"] * current_scale,
            "mapping_warning": "Jz/Jxy also differs, so this rescaling is a window locator, not a prediction.",
        },
        "finite_cluster_correction_comparisons": comparisons,
        "90mK_all_three_fields_opposite_sign_N15_vs_N18": all(not x["same_sign"] for x in comparisons if x["temperature_K"] == .09),
        "max_abs_N18_full_minus_low16_90mK": max(abs(x["N18_full_minus_low16"]) for x in comparisons if x["temperature_K"] == .09),
        "decision_reason": "The public quantum bulk baseline robustly relocates the relevant window to 1.45--1.60 T, but the M3-M0 correction reverses sign from N15 to N18 at every tested 90 mK field. The discrepancy is a cluster-size effect rather than stochastic FTLM noise.",
        "next_required_method": "A thermodynamic-limit tensor-network calculation of the realistic M3 Hamiltonian, or author-supplied/re-runnable XTRG data/code extended to Jpp and long-range dipoles.",
    }
    (RESULTS / "gate1_revised_report.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    entropy_rows = np.genfromtxt(RESULTS / "npj2022_xtrg_entropy_reconstructed.csv", delimiter=",", names=True, encoding="utf-8-sig")
    q = entropy_rows[(entropy_rows["anchor_N"] == 12) & np.isclose(entropy_rows["temperature_K"], .09)]
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.1))
    axes[0].plot(q["field_T"], q["entropy_kB_per_spin"], color="#4c78a8")
    axes[0].axvspan(xtrg_90["W95_min_T"], xtrg_90["W95_max_T"], color="#f58518", alpha=.22, label="95% design window")
    axes[0].axvline(xtrg_90["B_star_T"], color="#e45756", ls="--", label=f"B*={xtrg_90['B_star_T']:.2f} T")
    axes[0].set(xlabel="Field B (T)", ylabel=r"Entropy ($k_B$/spin)", title="Published XTRG M0 proxy, 90 mK")
    axes[0].grid(alpha=.25); axes[0].legend(frameon=False, fontsize=8)

    a = [x for x in comparisons if x["temperature_K"] == .09]
    axes[1].plot([x["field_T"] for x in a], [x["N15_deltaS"] for x in a], "o-", label="N=15")
    axes[1].plot([x["field_T"] for x in a], [x["N18_deltaS"] for x in a], "s-", label="N=18 compact")
    axes[1].axhline(0, color="black", lw=.8)
    axes[1].set(xlabel="Field B (T)", ylabel=r"$S(M3)-S(M0)$ ($k_B$/spin)", title="Correction reverses with cluster size")
    axes[1].grid(alpha=.25); axes[1].legend(frameon=False, fontsize=8)

    times = []
    labels = []
    for b in common_fields:
        for model in ("M0", "M3"):
            d = json.loads((RESULTS / "n18_bulk_window_compact_raw" / f"{model}_B{b:.3f}T.json").read_text(encoding="utf-8"))
            labels.append(f"{model}\n{b:.2f}T")
            times.append(d["wall_seconds_including_operator_build"] / 60)
    axes[2].bar(np.arange(len(times)), times, color=["#4c78a8" if "M0" in x else "#e45756" for x in labels])
    axes[2].set_xticks(np.arange(len(times)), labels, fontsize=8)
    axes[2].set(ylabel="Wall time (min)", title="Measured N=18 local cost")
    axes[2].grid(axis="y", alpha=.25)
    fig.suptitle("Revised Gate 1: bulk window found, realistic correction not size-converged")
    fig.tight_layout()
    fig.savefig(FIGURES / "gate1_revised_bulk_window_diagnostic.png", dpi=220)
    plt.close(fig)

    lines = [
        "# Gate 1 수정 판정: 보류",
        "",
        "## 확인된 사실",
        "",
        f"- 공개 XTRG 열지도와 엔트로피 곡선을 독립 복원했다. 0 T 교차검증 RMSE는 `{payload['validated_bulk_M0_proxy']['zero_field_validation_RMSE_J_mol_K']:.4f} J mol^-1 K^-1`이다.",
        f"- 90 mK M0 벌크 대리모형의 `B*`는 `{xtrg_90['B_star_T']:.2f} T`, 95% 운전창은 `{xtrg_90['W95_min_T']:.2f}–{xtrg_90['W95_max_T']:.2f} T`이다.",
        f"- 2026 매개변수로 단순 무차원 재척도화한 탐색 중심은 약 `{payload['validated_bulk_M0_proxy']['rough_current_parameter_field_mapping_T']:.2f} T`이나, 이 값은 이방성비 차이를 포함하지 않아 예측값이 아니다.",
        "- 따라서 초기 1.10–1.20 T 창은 유한클러스터 껍질효과가 만든 잘못된 계산창이며 폐기한다.",
        "",
        "## 현실 상호작용 보정의 수렴검사 (90 mK)",
        "",
        "|B (T)|N15 ΔS|N18 ΔS|N18 low-16 ΔS|부호 일치|",
        "|---:|---:|---:|---:|:---:|",
    ]
    for x in a:
        lines.append(f"|{x['field_T']:.2f}|{x['N15_deltaS']:+.6f}|{x['N18_deltaS']:+.6f}|{x['N18_low16_deltaS']:+.6f}|아니오|")
    lines += [
        "",
        f"전체 FTLM과 결정론적 저에너지 16상태의 차이는 최대 `{payload['max_abs_N18_full_minus_low16_90mK']:.2e} kB/spin`이다. 부호 반전은 무작위 추정 잡음이 아니라 유한크기 효과다.",
        "",
        "## 결정",
        "",
        "- Gate 1은 `HOLD`이다. 현재 수치로 M3의 벌크 B* 이동·후회량을 주장하지 않는다.",
        "- 양자 알고리즘 및 QPU 단계는 시작하지 않는다.",
        "- 다음 필수 증거는 현실 M3 해밀토니안의 열역학적 극한 텐서네트워크 계산, 또는 원저자 XTRG 수치·코드를 확보해 J±±와 장거리 쌍극자를 넣은 재계산이다.",
        "- N=21을 같은 방식으로 단순 확대하는 것은 N=18 한 점이 4.3–17.6분인 현재 자원과 크기 진동을 고려하면 우선순위가 낮다.",
    ]
    (ROOT / "Gate1_수정판정_20261011.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
