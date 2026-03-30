import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import build_matrix_promotion_family_status
from src.presentation import render_matrix_promotion_family_status_markdown


def test_matrix_promotion_readiness_marks_pea_and_soy_candidate_sets_ready_but_not_externally_unlocked(matrix_promotion_family_rows):
    rows = matrix_promotion_family_rows
    markdown = render_matrix_promotion_family_status_markdown(rows)
    by_protein = {row.protein_type: row for row in rows}

    assert by_protein["pea_iso"].candidate_set_ready is True
    assert by_protein["pea_iso"].external_assessment_unlocked is False
    assert by_protein["soy_iso"].candidate_set_ready is True
    assert by_protein["soy_iso"].external_assessment_unlocked is False
    assert "missing external meaty-positive benchmark" in by_protein["pea_iso"].blocker
    assert "Matrix Promotion Family Readiness" in markdown