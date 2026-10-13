from pathlib import Path

import numpy as np
import pytest

from engine import atomic_io


def _temp_files(directory: Path, prefix: str) -> list[Path]:
    return list(directory.glob(f".{prefix}.*"))


def test_atomic_text_write_replaces_complete_file_in_same_directory(tmp_path, monkeypatch):
    target = tmp_path / "campaign_summary.json"
    target.write_text('{"old": true}', encoding="utf-8")
    real_replace = atomic_io.os.replace
    seen = []

    def checked_replace(source, destination):
        source_path = Path(source)
        destination_path = Path(destination)
        assert source_path.parent == destination_path.parent
        seen.append((source_path, destination_path))
        real_replace(source, destination)

    monkeypatch.setattr(atomic_io.os, "replace", checked_replace)
    atomic_io.atomic_write_text(target, '{"new": true}')

    assert target.read_text(encoding="utf-8") == '{"new": true}'
    assert len(seen) == 1
    assert _temp_files(tmp_path, target.name) == []


def test_atomic_text_replace_failure_preserves_old_artifact_and_cleans_temp(
    tmp_path, monkeypatch
):
    target = tmp_path / "ledgers.json"
    target.write_text('{"complete": true}', encoding="utf-8")

    def fail_replace(source, destination):
        raise OSError("injected replace failure")

    monkeypatch.setattr(atomic_io.os, "replace", fail_replace)
    with pytest.raises(OSError, match="injected replace failure"):
        atomic_io.atomic_write_text(target, '{"partial":')

    assert target.read_text(encoding="utf-8") == '{"complete": true}'
    assert _temp_files(tmp_path, target.name) == []


def test_atomic_npz_checkpoint_round_trips_and_replaces_existing_file(
    tmp_path, monkeypatch
):
    target = tmp_path / "checkpoint.npz"
    atomic_io.atomic_savez_compressed(target, time=np.asarray(1.25), field=np.arange(4))
    with np.load(target) as saved:
        assert float(saved["time"]) == pytest.approx(1.25)
        np.testing.assert_array_equal(saved["field"], np.arange(4))

    real_replace = atomic_io.os.replace
    seen = []

    def checked_replace(source, destination):
        assert Path(source).parent == Path(destination).parent
        seen.append(Path(source))
        real_replace(source, destination)

    monkeypatch.setattr(atomic_io.os, "replace", checked_replace)
    atomic_io.atomic_savez_compressed(target, time=np.asarray(2.5), field=np.arange(3))

    with np.load(target) as saved:
        assert float(saved["time"]) == pytest.approx(2.5)
        np.testing.assert_array_equal(saved["field"], np.arange(3))
    assert len(seen) == 1
    assert _temp_files(tmp_path, target.stem) == []


def test_atomic_npz_replace_failure_preserves_previous_checkpoint(tmp_path, monkeypatch):
    target = tmp_path / "checkpoint.npz"
    atomic_io.atomic_savez_compressed(target, time=np.asarray(4.0), field=np.arange(2))

    def fail_replace(source, destination):
        raise OSError("injected checkpoint replace failure")

    monkeypatch.setattr(atomic_io.os, "replace", fail_replace)
    with pytest.raises(OSError, match="injected checkpoint replace failure"):
        atomic_io.atomic_savez_compressed(
            target, time=np.asarray(9.0), field=np.arange(5)
        )

    with np.load(target) as saved:
        assert float(saved["time"]) == pytest.approx(4.0)
        np.testing.assert_array_equal(saved["field"], np.arange(2))
    assert _temp_files(tmp_path, target.stem) == []


def test_atomic_npz_requires_npz_suffix(tmp_path):
    with pytest.raises(ValueError, match=r"end in .npz"):
        atomic_io.atomic_savez_compressed(
            tmp_path / "checkpoint.bin", value=np.arange(1)
        )
