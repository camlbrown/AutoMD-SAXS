"""Unit tests for ``automd_saxs.paths``. No external deps; runs standalone."""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs.paths import JobPaths  # noqa: E402


def test_simulation_dir_naming():
    p = JobPaths("/jobs", "1AKI")
    assert p.simulation_dir == os.path.join(os.path.abspath("/jobs"), "1AKI_simulation")


def test_stage_dirs_present():
    p = JobPaths("/jobs", "1AKI")
    assert p.pdb2gmx_dir.endswith("1AKI_simulation/pdb2gmx")
    assert p.genion_dir.endswith("1AKI_simulation/genion")
    assert p.saxs_dir.endswith("1AKI_simulation/SAXS")


def test_repeat_dirs_follow_n_repeats():
    p = JobPaths("/jobs", "x", n_repeats=2)
    assert p.repeat_dir(1).endswith("production/rep1")
    assert p.repeat_dir(2).endswith("production/rep2")
    try:
        p.repeat_dir(3)
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for out-of-range repeat")


def test_all_dirs_counts():
    # 1 sim dir + 10 stage dirs + 3 dirs per repeat
    p = JobPaths("/jobs", "x", n_repeats=3)
    assert len(p.all_dirs()) == 1 + 10 + 3 * 3
    # parents precede children
    dirs = p.all_dirs()
    assert dirs.index(p.production_dir) < dirs.index(p.repeat_dir(1))
    assert dirs.index(p.repeat_dir(1)) < dirs.index(p.extract_frames_dir(1))


def test_create_makes_dirs_and_is_idempotent():
    with tempfile.TemporaryDirectory() as tmp:
        p = JobPaths(tmp, "demo", n_repeats=1)
        made = p.create()
        for d in made:
            assert os.path.isdir(d)
        # second call must not raise
        p.create()


def test_rejects_zero_repeats():
    try:
        JobPaths("/jobs", "x", n_repeats=0)
    except ValueError:
        return
    raise AssertionError("expected ValueError")


def _run_standalone():
    funcs = sorted((n, o) for n, o in globals().items()
                   if n.startswith("test_") and callable(o))
    failures = []
    for name, fn in funcs:
        try:
            fn()
            print("PASS", name)
        except Exception as exc:  # noqa: BLE001
            failures.append((name, exc))
            print("FAIL", name, "->", exc)
    print("\n{0}/{1} passed".format(len(funcs) - len(failures), len(funcs)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_run_standalone())
