import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from scripts import import_react_ot_colab_artifacts as importer


def _build_archive(archive_path: Path) -> None:
    with ZipFile(archive_path, "w", compression=ZIP_DEFLATED) as archive:
        archive.writestr(
            "results/react_ot_pilot_manifest.json",
            json.dumps({"status": "completed", "targets": [{"target": "lysinoalanine_crosslink", "status": "ok"}]}),
        )
        archive.writestr(
            "results/lysinoalanine_crosslink_react_ot_seed.xyz",
            "2\nts\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n",
        )
        archive.writestr(
            "results/lysinoalanine_crosslink_react_ot_seed.json",
            json.dumps({"target": "lysinoalanine_crosslink", "status": "ok"}),
        )
        archive.writestr("results/ignored.txt", "ignore me")


def test_import_archive_extracts_only_repo_compatible_outputs(tmp_path, monkeypatch):
    monkeypatch.setattr(importer, "ROOT", tmp_path)
    archive_path = tmp_path / "react_ot_colab_artifacts.zip"
    out_dir = tmp_path / "results/computational_gap_refinement"
    _build_archive(archive_path)

    receipt = importer.import_archive(archive_path, out_dir)

    assert receipt["status"] == "imported"
    assert sorted(receipt["imported_files"]) == [
        "results/computational_gap_refinement/lysinoalanine_crosslink_react_ot_seed.json",
        "results/computational_gap_refinement/lysinoalanine_crosslink_react_ot_seed.xyz",
        "results/computational_gap_refinement/react_ot_pilot_manifest.json",
    ]
    assert receipt["ignored_files"] == ["results/ignored.txt"]
    assert (out_dir / "react_ot_pilot_manifest.json").exists()
    assert (out_dir / "lysinoalanine_crosslink_react_ot_seed.xyz").exists()
    receipt_path = out_dir / importer.RECEIPT_NAME
    assert receipt_path.exists()
    saved_receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    assert saved_receipt["status"] == "imported"


def test_import_archive_requires_manifest(tmp_path, monkeypatch):
    monkeypatch.setattr(importer, "ROOT", tmp_path)
    archive_path = tmp_path / "missing_manifest.zip"
    out_dir = tmp_path / "results/computational_gap_refinement"
    with ZipFile(archive_path, "w", compression=ZIP_DEFLATED) as archive:
        archive.writestr("results/lysinoalanine_crosslink_react_ot_seed.xyz", "2\nts\nH 0 0 0\nH 0 0 1\n")

    try:
        importer.import_archive(archive_path, out_dir)
    except ValueError as exc:
        assert "react_ot_pilot_manifest.json" in str(exc)
    else:
        raise AssertionError("expected missing-manifest archive to fail")