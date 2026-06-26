"""Subprocess wrapper with dry-run support.

External scientific binaries (GROMACS, ATSAS, Slurm) are treated as unavailable
by default. :class:`CommandRunner` lets the rest of the codebase be exercised
without them:

* ``dry_run=True`` records each planned command and runs nothing.
* ``dry_run=False`` checks the executable exists (raising
  :class:`MissingDependencyError` with an actionable message if not) before
  delegating to :mod:`subprocess`.

Commands are always passed as ``argv`` lists, never shell strings, so user
inputs are never interpolated into a shell.
"""

import shutil
from typing import List, Optional


class MissingDependencyError(RuntimeError):
    """Raised when a required external executable is not on ``PATH``."""


class PlannedCommand:
    """A recorded command (what a dry run would have executed)."""

    def __init__(self, label, argv, cwd=None, stdin=None):
        self.label = label
        self.argv = list(argv)
        self.cwd = cwd
        self.stdin = stdin

    def as_shell(self) -> str:
        """Best-effort human-readable rendering (display only, never executed)."""
        import shlex

        text = " ".join(shlex.quote(a) for a in self.argv)
        if self.stdin is not None:
            text = "echo {0} | {1}".format(shlex.quote(self.stdin), text)
        return text

    def __repr__(self):
        return "PlannedCommand(label={0!r}, argv={1})".format(self.label, self.argv)


class CommandRunner:
    """Plans or executes commands depending on ``dry_run``."""

    def __init__(self, dry_run=True):
        self.dry_run = dry_run
        self.planned = []  # type: List[PlannedCommand]

    def run(self, argv, label=None, cwd=None, stdin=None, check=True, capture=False):
        """Plan (dry run) or execute a command.

        Returns the :class:`PlannedCommand` in dry-run mode, or the
        ``subprocess.CompletedProcess`` otherwise. With ``capture=True`` the
        executed command's stdout is captured (used for ``sbatch --parsable``).
        """
        planned = PlannedCommand(label or (argv[0] if argv else ""), argv, cwd, stdin)
        self.planned.append(planned)
        if self.dry_run:
            return planned

        executable = argv[0]
        if shutil.which(executable) is None:
            raise MissingDependencyError(
                "Required executable {0!r} not found on PATH. Install it or run "
                "with dry_run=True to plan without executing.".format(executable)
            )

        import subprocess

        return subprocess.run(
            argv,
            cwd=cwd,
            input=(stdin.encode() if isinstance(stdin, str) else stdin),
            check=check,
            stdout=(subprocess.PIPE if capture else None),
        )

    def planned_labels(self) -> List[str]:
        return [p.label for p in self.planned]
