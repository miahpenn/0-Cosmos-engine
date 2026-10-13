"""Atomic file output for campaign evidence and restart checkpoints.

Writes are staged in the destination directory and made visible with os.replace,
so readers see either the previous complete artifact or the new complete one.
When supported, a best-effort parent-directory fsync is attempted after replacement;
filesystem and hardware behavior still limit power-loss durability guarantees.
This is output handling only; it does not change campaign or evolution semantics.
"""
import os
import tempfile
from pathlib import Path

import numpy as np


def atomic_write_text(path: str | Path, content: str, *, encoding: str = "utf-8") -> None:
    """Atomically replace a text artifact after flushing its complete contents."""
    target = Path(path)
    if not isinstance(content, str):
        raise TypeError("content must be a string")
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".tmp", dir=target.parent
    )
    temporary = Path(temp_name)
    try:
        with os.fdopen(fd, "w", encoding=encoding, newline="") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        _preserve_existing_mode(target, temporary)
        os.replace(temporary, target)
        _fsync_parent_directory(target.parent)
    except BaseException:
        try:
            temporary.unlink(missing_ok=True)
        except OSError:
            pass
        raise


def atomic_savez_compressed(path: str | Path, **arrays) -> None:
    """Atomically replace a NumPy compressed checkpoint."""
    target = Path(path)
    if target.suffix.lower() != ".npz":
        raise ValueError("compressed NumPy checkpoint path must end in .npz")
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(
        prefix=f".{target.stem}.", suffix=".npz", dir=target.parent
    )
    os.close(fd)
    temporary = Path(temp_name)
    try:
        np.savez_compressed(temporary, **arrays)
        with temporary.open("rb+") as stream:
            os.fsync(stream.fileno())
        _preserve_existing_mode(target, temporary)
        os.replace(temporary, target)
        _fsync_parent_directory(target.parent)
    except BaseException:
        try:
            temporary.unlink(missing_ok=True)
        except OSError:
            pass
        raise


def _preserve_existing_mode(target: Path, temporary: Path) -> None:
    """Keep an existing artifact's permissions when replacing it."""
    try:
        mode = target.stat().st_mode & 0o777
    except FileNotFoundError:
        return
    os.chmod(temporary, mode)


def _fsync_parent_directory(parent: Path) -> None:
    """Best-effort sync of rename metadata on platforms that support it."""
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    try:
        directory_fd = os.open(parent, flags)
    except OSError:
        return
    try:
        os.fsync(directory_fd)
    except OSError:
        # Some filesystems/platforms do not permit directory fsync.
        pass
    finally:
        try:
            os.close(directory_fd)
        except OSError:
            pass
