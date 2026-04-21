#!/usr/bin/env python3
# Designed for Linux

from __future__ import annotations

import os
from pathlib import Path
import subprocess
from typing import Any


ROOT = Path(__file__).resolve().parent
MYCOSV_BRIDGE_CPP = ROOT / "mycosv_cli_bridge.cpp"


def _resolve_executable(cmd: list[str]) -> list[str]:
    if not cmd:
        return cmd
    head = Path(cmd[0])
    if not head.exists() and head.suffix.lower() == ".exe":
        alt = head.with_suffix("")
        if alt.exists():
            fixed = cmd.copy()
            fixed[0] = str(alt)
            return fixed
    return cmd


def run_checked(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    cmd = _resolve_executable(cmd)
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        text=True,
        capture_output=True,
        check=True,
    )


def run_mycosv_command(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return run_checked(cmd, cwd=cwd)
