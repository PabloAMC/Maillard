from src.authority_benchmark_data import (
    load_irc_validation_cases,
    load_phase33_barrier_benchmarks,
    load_phase35_double_hybrid_benchmarks,
)


def test_load_phase33_barrier_benchmarks():
    benchmarks = load_phase33_barrier_benchmarks()

    assert len(benchmarks) == 9
    assert benchmarks["amadori_rearrangement"]["wb97mv_kcal_mol"] == 20.1
    assert benchmarks["trapping_hexanal"]["literature"]["best"] == 18.0


def test_load_phase35_double_hybrid_benchmarks():
    phase33 = load_phase33_barrier_benchmarks()
    phase35 = load_phase35_double_hybrid_benchmarks()

    assert set(phase35) == set(phase33)
    assert phase35["strecker_decarboxylation"]["revdsd_pbep86_d4_kcal_mol"] == 24.8
    assert phase35["amadori_rearrangement"]["wb97mv_kcal_mol"] == phase33["amadori_rearrangement"]["wb97mv_kcal_mol"]


def test_load_irc_validation_cases():
    cases = load_irc_validation_cases()
    case_ids = {case["case_id"] for case in cases}

    assert case_ids == {"amadori_proxy", "strecker_proxy"}
    assert all(case["ts_xyz"].startswith(("4\n", "5\n")) for case in cases)
    assert all(len(case["energies"]) == 3 for case in cases)