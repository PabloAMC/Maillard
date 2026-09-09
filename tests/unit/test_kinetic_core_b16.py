"""Wave B16 (2026-09-07): the thiol sinks against Schieberle 2000's 100 C series and Zhai 2021's TTCA decay."""
from __future__ import annotations

import pytest

from tests.support import wave_generator


def test_the_generator_splices_ten_rows_onto_b9_and_restores():
    with wave_generator("generate_kinetic_core_b16_fit") as B16:
        import generate_kinetic_core_b2_3_fit as B23

        B16.configure(False)
        ids = [r["id"] for r in B23.ACTIVE_FIT_ROWS]
        assert len(ids) == 64 and len(B16.B16_FIT_ROWS) == 10
        assert all(r["id"] in ids for r in B16.B16_FIT_ROWS)
        by = {r["id"]: r for r in B16.B16_FIT_ROWS}
        # Table IV ratios, verbatim
        assert by["schieberle_MFT_fold_720_over_30"]["target"] == pytest.approx(179.0 / 4.5)
        assert by["schieberle_FFT_fold_360_over_30"]["target"] == pytest.approx(110.0 / 2.0)
        assert by["schieberle_MFT_145C20min_over_100C360min"]["target"] == pytest.approx(1.0 / 13.0)
        assert by["schieberle_MFT_145C20min_over_100C360min"]["system"] == "hofmann_pentose_pH5"
        # the 100 C systems are the Hofmann pot, same charge and buffer, at 100 C
        s = B23.SYSTEMS["schieberle_100C_720"]
        assert s["initial"] == B23.SYSTEMS["hofmann_pentose_pH5"]["initial"] and s["t_c"] == 100.0 and s["minutes"] == 720.0
        assert s["buffer"] is B23.SYSTEMS["hofmann_pentose_pH5"]["buffer"]
        # Zhai 2021: c = c0 - k t at 60 min
        assert by["zhai2021_ttca_remaining_120C_60min"]["target"] == pytest.approx(9.9718 - 0.0651 * 60)
        assert B23.SYSTEMS["zhai2021_ttca_140"]["initial"]["TTCA"] == 10.0
        # no new coordinate; the ceiling is kept unless lifted
        lo, hi = B16.full_bounds()
        assert len(B16.FREE_KEYS) == 23 and hi[B16.THIOL_SINK_SLOT] == pytest.approx(102.0)
        B16.configure(True)
        lo, hi = B16.full_bounds()
        assert hi[B16.THIOL_SINK_SLOT] == pytest.approx(160.0)
        assert B16.OUT_FIT_REPORT.name == "kinetic_core_b16_lift_fit_report.json"
        B16.restore()
        assert len(B23.ACTIVE_FIT_ROWS) == 54 and "schieberle_100C_720" not in B23.SYSTEMS
