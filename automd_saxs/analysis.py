"""SAXS analysis orchestration.

Collects the outputs of the (ATSAS) SAXS-fitting stage from an analysis
directory, ranks structures by chi^2, writes a machine-readable summary table,
and folds the headline metrics into the run manifest.

Parsing lives in :mod:`automd_saxs.atsas`; this module is the orchestration layer
that turns parsed records into summary artefacts. It needs no ATSAS binaries --
it reads files that a previous (real) run produced, so it is fully testable from
fixtures. The same shape will host FoXS/MultiFoXS results in Phase 2.
"""

import csv
import glob
import os
from typing import Any, Dict, List, Optional

from . import atsas
from .manifest import Manifest, STATUS_COMPLETED


def collect_saxs_results(analysis_dir: str) -> Dict[str, Any]:
    """Gather parsed SAXS results from ``analysis_dir``.

    Looks for ``crysol_summary.txt`` (per-model chi^2/Rg), the first ``*.fit``
    (ensemble/best fit curve), and the first ``Rg_distr*.txt`` (ensemble Rg
    distribution). Missing files are simply absent from the result.
    """
    results = {
        "records": [],            # type: List[Dict[str, object]]
        "summary": atsas.summarize_fits([]),
        "rgDistribution": None,   # type: Optional[list]
        "fit": None,              # type: Optional[list]
        "filesFound": [],         # type: List[str]
    }  # type: Dict[str, Any]

    summary_path = os.path.join(analysis_dir, "crysol_summary.txt")
    if os.path.isfile(summary_path):
        with open(summary_path) as fh:
            records = atsas.parse_crysol_summary(fh.read())
        results["records"] = records
        results["summary"] = atsas.summarize_fits(records)
        results["filesFound"].append(summary_path)

    fits = sorted(glob.glob(os.path.join(analysis_dir, "*.fit")))
    if fits:
        with open(fits[0]) as fh:
            results["fit"] = atsas.parse_fit(fh.read())
        results["fitFile"] = fits[0]
        results["filesFound"].append(fits[0])

    rg_files = sorted(glob.glob(os.path.join(analysis_dir, "Rg_distr*.txt")))
    if rg_files:
        with open(rg_files[0]) as fh:
            results["rgDistribution"] = atsas.parse_rg_distribution(fh.read())
        results["rgDistFile"] = rg_files[0]
        results["filesFound"].append(rg_files[0])

    return results


def write_summary_csv(records: List[Dict[str, object]], path: str) -> str:
    """Write per-model ``frame,rg,chi2`` rows sorted by chi^2 (best first)."""
    ordered = sorted(records, key=lambda r: r["chi2"]) if records else []
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["rank", "frame", "model", "rg", "chi2"])
        for rank, rec in enumerate(ordered, start=1):
            writer.writerow([rank, rec["frame"], rec["model"], rec["rg"], rec["chi2"]])
    return path


def planned_commands(saxs_dat: Optional[str], smax: Optional[float] = None):
    """The ATSAS commands this analysis represents (for manifest traceability)."""
    if not saxs_dat:
        return []
    return [
        ("shanum", atsas.shanum_command(saxs_dat)),
        ("crysol", atsas.crysol_command(saxs_dat, smax=smax)),
        ("gajoe", atsas.gajoe_command("profiles.txt")),
    ]


def build_analysis_manifest(
    results: Dict[str, Any],
    config=None,
    summary_csv: Optional[str] = None,
    saxs_dat: Optional[str] = None,
) -> Manifest:
    """Build a ``completed`` manifest from collected SAXS results."""
    manifest = Manifest(config, status=STATUS_COMPLETED)
    summary = results.get("summary") or {}
    for key in ("bestChi2", "bestFrame", "rgMean"):
        if key in summary:
            manifest.set_metric(key, summary[key])
    manifest.set_parameter("nModelsScored", len(results.get("records", [])))

    if summary_csv:
        manifest.add_output("summaryTables", summary_csv)
    if results.get("fitFile"):
        manifest.add_output("saxsFits", results["fitFile"])
    if results.get("rgDistFile"):
        manifest.add_output("summaryTables", results["rgDistFile"])

    for name, command in planned_commands(saxs_dat):
        manifest.add_step(name, command=command)

    if not results.get("records"):
        manifest.add_note("no crysol_summary.txt found; metrics are empty")
    return manifest
