"""Executable documentation macros for the Zensical site."""

from __future__ import annotations

import os
import re
import subprocess
import tempfile
from pathlib import Path

from zensical.extensions.macros import MacroEnv

_PROJECT_ROOT = Path(__file__).parent
_WORK_DIRECTORY = tempfile.TemporaryDirectory(prefix="protein-fasta-cli-docs-")
_WORK_PATH = Path(_WORK_DIRECTORY.name)
_ANSI_ESCAPE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")
_LOGURU_PREFIX = re.compile(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3} \| [A-Z]+\s+\| .*? - ")


def define_env(environment: MacroEnv) -> None:
    """Register the shell-command macro used by the CLI walkthrough."""
    environment.macro(cli)


def cli(command: str) -> str:
    """Execute one authored Bash command and return a stable console transcript."""
    command = command.strip()
    process_environment = os.environ.copy()
    process_environment.update({"NO_COLOR": "1", "WORK": str(_WORK_PATH)})
    completed = subprocess.run(
        ["/bin/bash", "--noprofile", "--norc", "-euo", "pipefail", "-c", command],
        cwd=_PROJECT_ROOT,
        env=process_environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    output = _stable_output(completed.stdout)
    if completed.returncode:
        raise RuntimeError(
            f"CLI documentation command failed with exit code {completed.returncode}:\n"
            f"{command}\n{output}"
        )
    transcript = f"$ {command}"
    if output:
        transcript = f"{transcript}\n{output}"
    return f"```console\n{transcript}\n```"


def _stable_output(output: str) -> str:
    """Remove volatile timestamps, ANSI escapes, and temporary absolute paths."""
    lines: list[str] = []
    for raw_line in output.splitlines():
        line = _ANSI_ESCAPE.sub("", raw_line)
        line = _LOGURU_PREFIX.sub("", line)
        lines.append(line.replace(str(_WORK_PATH), "$WORK"))
    return "\n".join(lines).strip()
