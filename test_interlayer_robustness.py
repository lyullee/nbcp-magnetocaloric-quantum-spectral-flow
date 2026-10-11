import json
from pathlib import Path

import numpy as np

from bulk_polarized_stability import (
    EwaldSettings,
    InterlayerParameters,
    locate_critical_mode,
)
from bulk_to_engineering_surrogate import transform_field
from compile_refrigerator_design_chart import CELL_VOLUME_A3, DENSITY_G_CM3
from screen_interlayer_parameter_ensemble import (
    canonical_soft_mode,
    out_of_plane_bandwidth_microev,
)
from stress_test_thermodynamic_shape import deformed_profile


def run():
    m0 = canonical_soft_mode([1 / 3, 2 / 3, 0.0])
    assert m0["delta_inplane_rlu"] < 1e-12

    representative = InterlayerParameters()
    m4 = locate_critical_mode("M4", EwaldSettings(),
                              interlayer=representative)
    mode = canonical_soft_mode(m4["critical_wavevector_rlu"])
    assert abs(mode["delta_inplane_rlu"] - 0.02482641) < 2e-7
    assert abs(mode["qz_folded_rlu"] - 0.33428861) < 2e-7
    assert abs(out_of_plane_bandwidth_microev(representative) - 4.4) < 0.03

    fields = np.array([1.3, 1.4, 1.5])
    heat = np.array([0.2, 0.5, 0.3])
    source_hc, target_hc = 1.75, 1.65
    warped, values = deformed_profile(
        fields, heat, source_hc, target_hc, "critical_offset",
        peak_shift_t=0.0, left_scale=1.0, right_scale=1.0,
        heat_scale=1.0)
    assert np.allclose(warped, transform_field(
        fields, source_hc, target_hc, "critical_offset"))
    assert np.array_equal(values, heat)

    ion_density = 1.0 / (CELL_VOLUME_A3 * 1e-24)
    assert abs(ion_density / 5.8e21 - 1.0) < 0.02
    assert 4.1 < DENSITY_G_CM3 < 4.3

    root = Path(__file__).resolve().parent
    measurement = json.loads(
        (root / "results" / "peak_measurement_target_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert measurement["largest_symmetric_peak_uncertainty_with_common_W95_T"] == 0.002
    assert measurement["recommended_no_fit_sweep"]["point_count"] == 26
    assert measurement["recommended_no_fit_sweep"]["step_T"] == 0.004
    print("interlayer robustness validation passed")


if __name__ == "__main__":
    run()
