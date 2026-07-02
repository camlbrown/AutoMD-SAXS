"""Structural ensemble clustering / PCA orchestration.

Ports the legacy ``run_structural.py`` + ``structural_utils.py`` (CLoNe driver,
cluster summaries) into the package. Split into two layers:

* **Pure-stdlib helpers** -- cluster summary statistics (per-dimension centre,
  median, IQR; cluster sizes), CSV writing, and parsers for the legacy
  ``PCA_coords.txt`` / ``results.txt`` outputs. These are tested without numpy.
* **Numerical orchestration** -- :func:`run_clustering` loads a trajectory
  (mdtraj), optionally reduces with PCA (scikit-learn), runs :class:`CLoNe`, and
  writes a summary. Heavy deps are imported lazily with a clear error.

The clustering maths lives in :mod:`automd_saxs.clone` and is unchanged from the
legacy implementation, preserving the science.
"""

import csv
import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

from .command_runner import MissingDependencyError


# --------------------------------------------------------------------------- #
# Pure-stdlib statistics
# --------------------------------------------------------------------------- #

def _percentile(values: Sequence[float], q: float) -> float:
    """Linear-interpolation percentile matching numpy's default ``np.percentile``."""
    xs = sorted(values)
    n = len(xs)
    if n == 0:
        raise ValueError("percentile of empty sequence")
    if n == 1:
        return float(xs[0])
    rank = (q / 100.0) * (n - 1)
    lo = int(math.floor(rank))
    hi = int(math.ceil(rank))
    if lo == hi:
        return float(xs[lo])
    return float(xs[lo] + (xs[hi] - xs[lo]) * (rank - lo))


def cluster_summary(data, labels, centers, headers=None, labels_all=None) -> List[Dict]:
    """Per-cluster, per-dimension summary (centre value, median, IQR, sizes).

    Port of ``structural_utils.show_cluster_info`` returning structured data
    instead of printing a table. ``data`` is an indexable of rows (list of lists
    or a numpy array). ``labels`` are post-outlier cluster ids (-1 = outlier);
    ``centers`` are data indices, one per cluster. ``labels_all`` (optional,
    pre-outlier) is used for the element count, matching the legacy distinction.
    """
    n_dims = len(data[0]) if len(data) else 0
    if headers is None:
        headers = ["C{0}".format(d) for d in range(n_dims)]
    labels = list(labels)
    out = []
    for cl in range(len(centers)):
        members = [i for i, lab in enumerate(labels) if lab == cl]
        n_after = len(members)
        if labels_all is not None:
            n_all = sum(1 for lab in labels_all if lab == cl)
        else:
            n_all = n_after
        ci = centers[cl]
        dims = []
        for d in range(n_dims):
            center_val = float(data[ci][d])
            col = [float(data[i][d]) for i in members]
            if col:
                q1 = _percentile(col, 25)
                med = _percentile(col, 50)
                q3 = _percentile(col, 75)
                iqr = q3 - q1
            else:
                med = iqr = float("nan")
            dims.append({"header": headers[d], "center": center_val,
                         "median": med, "iqr": iqr})
        out.append({
            "cluster": cl,
            "center_index": int(ci),
            "n_elements": n_all,
            "n_after_outliers": n_after,
            "dimensions": dims,
        })
    return out


def write_cluster_summary_csv(summary: List[Dict], path: str) -> str:
    """Write the cluster summary in long form (one row per cluster/dimension)."""
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["cluster", "center_index", "n_elements", "n_after_outliers",
                         "dimension", "center", "median", "iqr"])
        for c in summary:
            for dim in c["dimensions"]:
                writer.writerow([c["cluster"] + 1, c["center_index"], c["n_elements"],
                                 c["n_after_outliers"], dim["header"],
                                 dim["center"], dim["median"], dim["iqr"]])
    return path


# --------------------------------------------------------------------------- #
# Parsers for legacy outputs
# --------------------------------------------------------------------------- #

def parse_pca_coords(text: str) -> Tuple[List[str], List[Tuple[float, ...]]]:
    """Parse ``PCA_coords.txt``: header ``PC1(0.43) PC2(0.21)`` then numeric rows."""
    lines = text.splitlines()
    if not lines:
        return [], []
    headers = lines[0].split()
    rows = []
    for line in lines[1:]:
        parts = line.split()
        if not parts:
            continue
        try:
            rows.append(tuple(float(p) for p in parts))
        except ValueError:
            continue
    return headers, rows


def parse_cluster_results(text: str) -> Dict[str, object]:
    """Parse the legacy ``results.txt``: centre indices line, then ``label core rho`` rows."""
    lines = [l for l in text.splitlines() if l.strip()]
    if not lines:
        return {"centers": [], "points": []}
    centers = [int(x) for x in lines[0].split()]
    points = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < 3:
            continue
        points.append((int(parts[0]), int(parts[1]), float(parts[2])))
    return {"centers": centers, "points": points}


# --------------------------------------------------------------------------- #
# Numerical orchestration (lazy deps)
# --------------------------------------------------------------------------- #

def run_clustering(
    traj: str,
    topo: str,
    out_dir: str,
    at_sel: str = "name CA",
    pca: int = 0,
    pdc: float = 4.0,
    n_resize: float = 4.0,
    filt: float = 0.1,
    verbose: bool = False,
) -> Dict[str, object]:
    """Load a trajectory, optionally PCA-reduce, run CLoNe, write a summary.

    Non-interactive port of ``run_structural.py``. Requires mdtraj + scikit-learn
    + numpy + scipy; raises :class:`MissingDependencyError` if any is absent.
    Returns a dict with the summary, labels, centres, and written output files.
    """
    try:
        import numpy as np
        import mdtraj
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
    except ImportError as exc:
        raise MissingDependencyError(
            "Structural clustering requires mdtraj, scikit-learn, numpy and scipy "
            "(all in automdsaxs.yml). Original import error: {0}".format(exc)
        )
    from .clone import CLoNe

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    struct_ens = mdtraj.load(traj, top=topo)

    # Clustering needs enough frames to be meaningful (and PCA needs at least
    # `pca`+1 samples). Skip cleanly with a clear reason for very short runs
    # rather than letting sklearn raise a cryptic error.
    min_frames = max(3, (pca + 1) if pca else 3)
    if struct_ens.n_frames < min_frames:
        return {
            "summary": [],
            "labels": [],
            "centers": [],
            "nClusters": 0,
            "outputFiles": [],
            "skipped": "too few frames for clustering ({0} < {1})".format(
                struct_ens.n_frames, min_frames),
        }

    selection = struct_ens.topology.select(at_sel)
    coords = struct_ens.xyz[:, selection]
    coords = coords.reshape(coords.shape[0], coords.shape[1] * coords.shape[2])
    headers = ["C{0}".format(x) for x in range(coords.shape[1])]
    original_coords = coords.copy()

    if pca:
        pca_obj = PCA(n_components=pca)
        coords = pca_obj.fit_transform(coords)
        ev = pca_obj.explained_variance_ratio_
        headers = ["PC{0}({1:.2f})".format(x + 1, ev[x]) for x in range(pca)]
        original_coords = coords
        with open(os.path.join(out_dir, "PCA_coords.txt"), "w") as fh:
            fh.write(" ".join(headers) + "\n")
            for row in coords:
                fh.write(" ".join("{0:f}".format(v) for v in row) + "\n")
    else:
        coords = StandardScaler().fit_transform(coords)

    clone = CLoNe(pdc=pdc, n_resize=n_resize, filt=filt, verbose=verbose)
    clone.fit(coords)

    summary = cluster_summary(
        [list(r) for r in np.asarray(original_coords)],
        list(clone.labels_), list(clone.centers), headers, list(clone.labels_all))
    csv_path = write_cluster_summary_csv(summary, os.path.join(out_dir, "Summary_clusters.csv"))

    return {
        "summary": summary,
        "labels": list(clone.labels_),
        "centers": list(clone.centers),
        "nClusters": len(clone.centers),
        "outputFiles": [csv_path],
    }


def run_clustering_sweep(
    traj: str,
    topo: str,
    out_dir: str,
    at_sel: str = "name CA",
    pca: int = 2,
    pdc_values: Sequence[float] = (1, 2, 3, 4, 5, 6, 7),
    n_resize: float = 4.0,
    filt: float = 0.1,
    verbose: bool = False,
) -> Dict[str, object]:
    """Run CLoNe at several ``pdc`` values over one PCA-reduced trajectory.

    The original AutoMD-SAXS ``run_CLoNe.sh`` clustered the aligned trajectory at
    a range of ``pdc`` values (1, 3, 5, 7) so the user could pick the granularity.
    PCA is independent of ``pdc``, so we compute it once (writing ``PCA_coords.txt``)
    and only re-run the cheap CLoNe step per ``pdc``. Returns the shared PCA coords
    plus, for each ``pdc``, the cluster labels / centre count / summary.

    Returns ``{"pcaCoords", "headers", "sweep": [{pdc, nClusters, labels,
    summary, outputFiles}], "skipped"?}``.
    """
    try:
        import numpy as np
        import mdtraj
        from sklearn.decomposition import PCA
    except ImportError as exc:
        raise MissingDependencyError(
            "Structural clustering requires mdtraj, scikit-learn, numpy and scipy "
            "(all in automdsaxs.yml). Original import error: {0}".format(exc)
        )
    from .clone import CLoNe

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    struct_ens = mdtraj.load(traj, top=topo)
    # PCA needs at least pca+1 samples; CLoNe needs a handful of frames.
    min_frames = max(3, pca + 1)
    if struct_ens.n_frames < min_frames:
        return {"pcaCoords": [], "headers": [], "sweep": [],
                "skipped": "too few frames for clustering ({0} < {1})".format(
                    struct_ens.n_frames, min_frames)}

    selection = struct_ens.topology.select(at_sel)
    coords = struct_ens.xyz[:, selection]
    coords = coords.reshape(coords.shape[0], coords.shape[1] * coords.shape[2])

    pca_obj = PCA(n_components=pca)
    reduced = pca_obj.fit_transform(coords)
    ev = pca_obj.explained_variance_ratio_
    headers = ["PC{0}({1:.2f})".format(x + 1, ev[x]) for x in range(pca)]
    with open(os.path.join(out_dir, "PCA_coords.txt"), "w") as fh:
        fh.write(" ".join(headers) + "\n")
        for row in reduced:
            fh.write(" ".join("{0:f}".format(v) for v in row) + "\n")

    reduced_rows = [list(r) for r in np.asarray(reduced)]
    sweep = []
    for pdc in pdc_values:
        try:
            clone = CLoNe(pdc=pdc, n_resize=n_resize, filt=filt, verbose=verbose)
            clone.fit(reduced)
            summary = cluster_summary(
                reduced_rows, list(clone.labels_), list(clone.centers),
                headers, list(clone.labels_all))
            csv_path = write_cluster_summary_csv(
                summary, os.path.join(out_dir, "Summary_clusters_pdc{0}.csv".format(pdc)))
            sweep.append({
                "pdc": pdc,
                "nClusters": len(clone.centers),
                "labels": [int(x) for x in clone.labels_],
                "summary": summary,
                "outputFiles": [csv_path],
            })
        except Exception as exc:  # noqa: BLE001 - one pdc failing must not kill the sweep
            sweep.append({"pdc": pdc, "nClusters": 0, "labels": [],
                          "summary": [], "outputFiles": [], "error": str(exc)})

    return {
        "pcaCoords": [[float(v) for v in row] for row in reduced_rows],
        "headers": headers,
        "sweep": sweep,
    }
