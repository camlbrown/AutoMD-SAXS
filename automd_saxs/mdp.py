"""GROMACS ``.mdp`` templating.

The legacy ``run_MD.sh`` edited the production ``.mdp`` *in place* with
``sed -i`` -- mutating a git-tracked template on every run and corrupting the
repository state (a non-idempotent side effect, and one of the original
"outright bugs"). This module instead renders a template **into a string / the
job directory**, leaving the tracked template untouched.

``render_mdp`` is pure (text in, text out) and preserves comments, ordering, and
alignment, replacing only the values of the keys you override. GROMACS treats
``-`` and ``_`` interchangeably and is case-insensitive in mdp option names, so
key matching here normalises both.
"""

import re
from typing import Dict, Optional

from .config import JobConfig

# key = value ; optional comment   (whitespace, tabs, and multi-token values ok)
_ASSIGN = re.compile(
    r"^(?P<lead>\s*)"
    r"(?P<key>[A-Za-z0-9_\-]+)"
    r"(?P<pre>\s*)=(?P<post>\s*)"
    r"(?P<val>.*?)"
    r"(?P<trail>\s*)"
    r"(?P<comment>;.*)?$"
)


def normalize_key(key: str) -> str:
    """Canonical form for mdp option names (case-insensitive, ``-``/``_`` merged)."""
    return key.strip().lower().replace("-", "_")


def parse_mdp(text: str) -> Dict[str, str]:
    """Parse an mdp into ``{normalized_key: value}`` (values stripped)."""
    result = {}  # type: Dict[str, str]
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue
        match = _ASSIGN.match(line)
        if not match or "=" not in line:
            continue
        result[normalize_key(match.group("key"))] = match.group("val").strip()
    return result


def render_mdp(text: str, overrides: Dict[str, object], append_missing: bool = True) -> str:
    """Return ``text`` with the values of ``overrides`` applied.

    Keys are matched case-insensitively with ``-``/``_`` equivalence. Comments,
    ordering, and surrounding whitespace are preserved. Override keys that are
    not present are appended at the end (when ``append_missing``) so the result
    always reflects every requested value.
    """
    wanted = {normalize_key(k): v for k, v in overrides.items()}
    applied = set()
    out_lines = []

    for line in text.splitlines():
        match = _ASSIGN.match(line)
        if match and "=" in line and not line.strip().startswith(";"):
            norm = normalize_key(match.group("key"))
            if norm in wanted:
                new_val = _fmt(wanted[norm])
                comment = match.group("comment") or ""
                rebuilt = "{lead}{key}{pre}={post}{val}".format(
                    lead=match.group("lead"),
                    key=match.group("key"),
                    pre=match.group("pre"),
                    post=match.group("post"),
                    val=new_val,
                )
                if comment:
                    rebuilt = rebuilt + " " + comment
                out_lines.append(rebuilt)
                applied.add(norm)
                continue
        out_lines.append(line)

    if append_missing:
        for norm, value in wanted.items():
            if norm not in applied:
                out_lines.append("{0} = {1}".format(norm, _fmt(value)))

    result = "\n".join(out_lines)
    if text.endswith("\n"):
        result += "\n"
    return result


def production_overrides(config: JobConfig) -> Dict[str, object]:
    """Overrides for the production mdp derived from the config.

    Generalises the legacy ``sed`` edit (which only touched ``nsteps``):

    * ``nsteps`` from :meth:`JobConfig.number_of_steps`;
    * ``dt`` (ps) from the timestep, kept consistent with ``nsteps``.
    """
    return {
        "nsteps": config.number_of_steps(),
        "dt": config.timestep_fs / 1000.0,
    }


def render_production_mdp(config: JobConfig, template_text: str) -> str:
    """Render the production mdp text for ``config`` (pure; no file I/O)."""
    return render_mdp(template_text, production_overrides(config))


def write_rendered_mdp(template_path: str, out_path: str, overrides: Dict[str, object]) -> str:
    """Read a template, render with ``overrides``, write to ``out_path``.

    The thin I/O wrapper around :func:`render_mdp`; the tracked template is read
    only, never modified.
    """
    with open(template_path, "r") as handle:
        text = handle.read()
    rendered = render_mdp(text, overrides)
    with open(out_path, "w") as handle:
        handle.write(rendered)
    return out_path


def _fmt(value: object) -> str:
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return ("%g" % value)
    return str(value)
