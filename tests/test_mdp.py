"""Unit tests for ``automd_saxs.mdp``. No external deps; runs standalone."""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automd_saxs import mdp  # noqa: E402
from automd_saxs.config import JobConfig  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD_TEMPLATE = os.path.join(REPO, "mdp_files", "md.mdp")

SAMPLE = (
    "; production\n"
    "integrator              = md        ; leap-frog\n"
    "nsteps                  = 500000 ; comment kept\n"
    "dt                      = 0.002     ; 2 fs\n"
    "tc-grps                 = Protein Non-Protein   ; two groups\n"
)


def test_parse_mdp_reads_values():
    values = mdp.parse_mdp(SAMPLE)
    assert values["nsteps"] == "500000"
    assert values["dt"] == "0.002"
    # hyphen/underscore equivalence
    assert values["tc_grps"] == "Protein Non-Protein"


def test_render_replaces_value_and_keeps_comment():
    out = mdp.render_mdp(SAMPLE, {"nsteps": 25000000})
    assert "comment kept" in out
    assert mdp.parse_mdp(out)["nsteps"] == "25000000"
    # other keys untouched
    assert mdp.parse_mdp(out)["dt"] == "0.002"


def test_render_key_normalisation():
    # override given with different case / separator still matches
    out = mdp.render_mdp(SAMPLE, {"NSTEPS": 7})
    assert mdp.parse_mdp(out)["nsteps"] == "7"


def test_render_appends_missing_key():
    out = mdp.render_mdp(SAMPLE, {"gen_seed": 12345})
    assert mdp.parse_mdp(out)["gen_seed"] == "12345"


def test_render_preserves_multitoken_values():
    out = mdp.render_mdp(SAMPLE, {"nsteps": 10})
    assert mdp.parse_mdp(out)["tc_grps"] == "Protein Non-Protein"


def test_production_overrides_from_config():
    cfg = JobConfig(protein_file="p.pdb", simulation_time_ns=50, timestep_fs=2.0)
    ov = mdp.production_overrides(cfg)
    assert ov["nsteps"] == 50 * 500000
    assert abs(ov["dt"] - 0.002) < 1e-12


def test_render_production_on_real_template():
    if not os.path.isfile(MD_TEMPLATE):
        print("skip: md.mdp template not found")
        return
    with open(MD_TEMPLATE) as fh:
        template = fh.read()
    cfg = JobConfig(protein_file="p.pdb", simulation_time_ns=20, timestep_fs=2.0)
    rendered = mdp.render_production_mdp(cfg, template)
    assert mdp.parse_mdp(rendered)["nsteps"] == "10000000"
    # template string itself is unchanged (render is pure)
    assert mdp.parse_mdp(template)["nsteps"] == "500000"


def test_write_rendered_does_not_touch_template():
    with tempfile.TemporaryDirectory() as tmp:
        tmpl = os.path.join(tmp, "md.mdp")
        with open(tmpl, "w") as fh:
            fh.write(SAMPLE)
        out = os.path.join(tmp, "rendered.mdp")
        mdp.write_rendered_mdp(tmpl, out, {"nsteps": 999})
        # output has the new value; template still has the old one
        assert mdp.parse_mdp(open(out).read())["nsteps"] == "999"
        assert mdp.parse_mdp(open(tmpl).read())["nsteps"] == "500000"


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
