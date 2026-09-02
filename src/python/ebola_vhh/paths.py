"""Shared path loader for a later Python rewrite of the figure code."""

from __future__ import annotations

import json
from pathlib import Path


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def load_paths() -> dict:
    root = repo_root()
    cfg = root / "config" / "paths.json"
    example = root / "config" / "paths.example.json"
    if not cfg.exists():
        cfg.write_text(example.read_text(encoding="utf-8"), encoding="utf-8")
    data = json.loads(cfg.read_text(encoding="utf-8"))
    data["repo_root"] = str(root)
    data["output"] = str(root / "output")
    data["manifest"] = str(root / "config" / "figure_manifest.json")
    return data


def load_manifest() -> dict:
    path = Path(load_paths()["manifest"])
    return json.loads(path.read_text(encoding="utf-8"))
