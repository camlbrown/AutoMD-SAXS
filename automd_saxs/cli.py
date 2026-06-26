"""Command-line entry points.

Phase 1 implements the dry-run planning contract from the project brief::

    python -m automd_saxs plan --config job.json [--pdb model.pdb] [--dry-run]

``plan`` validates the configuration and prints the planned scientific steps and
commands (GROMACS setup pipeline + Slurm submission DAG + job directory layout)
**without** running GROMACS, ATSAS, or Slurm. ``--config`` accepts either the new
``job.json`` or a legacy ``configurations.txt`` (auto-detected by extension).
"""

import argparse
import os
import sys

from . import analysis as analysis_mod
from . import dmax as dmax_mod
from . import gromacs_commands as gmx
from . import ligand as ligand_mod
from . import mdp as mdp_mod
from . import manifest as manifest_mod
from . import slurm as slurm_mod
from .command_runner import CommandRunner, MissingDependencyError
from .config import BoxShape, JobConfig, SystemType
from .legacy import config_from_legacy_file


def load_config(path: str) -> JobConfig:
    """Load a :class:`JobConfig` from ``job.json`` or legacy ``configurations.txt``."""
    if path.lower().endswith(".json"):
        with open(path, "r") as handle:
            return JobConfig.from_json(handle.read())
    return config_from_legacy_file(path)


def _protein_name(protein_file: str) -> str:
    base = os.path.basename(protein_file)
    return base[:-4] if base.lower().endswith(".pdb") else base


def cmd_plan(args) -> int:
    config = load_config(args.config)
    out = []

    out.append("AutoMD-SAXS plan (dry run)")
    out.append("=" * 60)
    out.append("Configuration")
    out.append("-" * 60)
    for line in config.to_json().splitlines():
        out.append("  " + line)

    out.append("")
    out.append("Effective parameters")
    out.append("-" * 60)
    out.append("  system            : {0}".format(config.system.value))
    out.append("  force field       : {0}".format(config.force_field_name))
    out.append("  water model       : {0}".format(config.water_model()))
    out.append("  production nsteps  : {0} ({1} ns at {2} fs)".format(
        config.number_of_steps(), config.simulation_time_ns, config.timestep_fs))
    out.append("  genion -conc       : {0} M".format(config.genion_concentration()))
    out.append("  repeats           : {0}".format(config.n_repeats))
    out.append("  uses SAXS         : {0}".format(config.uses_saxs))

    # --- Dmax / box padding (needs the structure) ------------------------
    box_padding = None
    out.append("")
    out.append("Box padding")
    out.append("-" * 60)
    if args.pdb and os.path.isfile(args.pdb):
        model = dmax_mod.model_dmax_nm(args.pdb)
        effective = dmax_mod.resolve_dmax_nm(config.dmax_nm, model)
        box_padding = dmax_mod.box_padding_nm(effective, model)
        out.append("  model Dmax        : {0:.3f} nm".format(model))
        out.append("  effective Dmax    : {0:.3f} nm".format(effective))
        out.append("  box padding (-d)  : {0:.3f} nm".format(box_padding))
    else:
        out.append("  (provide --pdb to compute model Dmax and box padding)")

    # --- GROMACS setup pipeline -----------------------------------------
    out.append("")
    out.append("GROMACS setup commands")
    out.append("-" * 60)
    setup_steps = []
    if config.system is SystemType.PROTEIN_LIGAND:
        ligand_names = []
        if args.pdb and os.path.isfile(args.pdb):
            with open(args.pdb) as fh:
                _, ligands = ligand_mod.split_complex(fh.read())
            ligand_names = ["ligand_{0}".format(num) for _, num, _ in ligands]
            out.append("  detected ligands  : {0}".format(
                ", ".join(ligand_names) or "(none found)"))
        else:
            ligand_names = ["ligand_<N>"]
            out.append("  (provide --pdb to detect ligands; showing a template)")
        for label, argv, stdout in ligand_mod.plan_ligand_setup(ligand_names, config.disulfide):
            redirect = "  > {0}".format(stdout) if stdout else ""
            out.append("  {0:16} {1}{2}".format(label, " ".join(argv), redirect))
        out.append("  (then box/solvate/genion proceed as the protein path, on the complex)")
        setup_steps = [(label, argv) for label, argv, _ in
                       ligand_mod.plan_ligand_setup(ligand_names, config.disulfide)]
    elif config.box_shape is BoxShape.AUTO:
        out.append("  [note] box_shape=auto resolves via a runtime gyrate "
                   "heuristic; cannot be planned statically.")
    else:
        ions_mdp = os.path.join(args.mdp_dir, "ions.mdp")
        if box_padding is None:
            # All steps except the box step, which needs the structure.
            setup_steps = [
                ("pdb2gmx", gmx.pdb2gmx_command(config)),
                ("editconf_box", [gmx.DEFAULT_GMX, "editconf", "-f", "GMX.gro",
                                  "-o", "1.gro", "-bt",
                                  config.box_shape.value, "-d", "<needs --pdb>"]),
                ("editconf_center", gmx.editconf_center_command()),
                ("solvate", gmx.solvate_command()),
                ("ions_grompp", gmx.ions_grompp_command(ions_mdp)),
                ("genion", gmx.genion_command(config)),
            ]
        else:
            setup_steps = gmx.setup_pipeline(config, box_padding, ions_mdp)
        for label, argv in setup_steps:
            out.append("  {0:16} {1}".format(label, " ".join(argv)))
        out.append("  (genion reads group {0!r} on stdin)".format(gmx.GENION_REPLACE_GROUP))

    # --- production mdp (rendered, never sed-edited in place) -------------
    out.append("")
    out.append("Production MDP (rendered into job dir, template left untouched)")
    out.append("-" * 60)
    overrides = mdp_mod.production_overrides(config)
    template_name = "lig_md.mdp" if config.system is SystemType.PROTEIN_LIGAND else "md.mdp"
    template_path = os.path.join(args.mdp_dir, template_name)
    out.append("  template          : {0}".format(template_path))
    out.append("  overrides         : nsteps={0}, dt={1} ps".format(
        overrides["nsteps"], overrides["dt"]))
    if os.path.isfile(template_path):
        with open(template_path) as fh:
            before = mdp_mod.parse_mdp(fh.read()).get("nsteps")
        out.append("  template nsteps    : {0} (unchanged on disk)".format(before))

    # --- Slurm submission DAG -------------------------------------------
    out.append("")
    out.append("Slurm submission DAG")
    out.append("-" * 60)
    paths = _build_paths(config, args)
    submissions = slurm_mod.plan_submissions(
        config, args.slurm_dir, paths.simulation_dir, config.partition_flag()
    )
    for job, argv in submissions:
        deps = ",".join(job.depends_on) if job.depends_on else "-"
        out.append("  {0:7} (after {1:12}) {2}".format(job.name, deps, " ".join(argv)))

    # --- directory layout ------------------------------------------------
    out.append("")
    out.append("Job directories ({0} would be created)".format(len(paths.all_dirs())))
    out.append("-" * 60)
    for d in paths.all_dirs():
        out.append("  " + d)

    # --- manifest --------------------------------------------------------
    manifest = manifest_mod.build_plan_manifest(config, box_padding, setup_steps)
    if args.manifest:
        manifest.write(args.manifest)
        out.append("")
        out.append("Manifest written to: {0}".format(args.manifest))

    print("\n".join(out))
    return 0


def _build_paths(config, args):
    from .paths import JobPaths

    base = args.work_dir or os.getcwd()
    return JobPaths(base, _protein_name(config.protein_file), config.n_repeats)


def cmd_dmax(args) -> int:
    """Print model/effective Dmax and box padding as machine-readable KEY=VALUE lines.

    Replaces the inline ``python3``/``bc`` block in ``run_MD.sh``; the shell can
    ``eval`` or grep the output.
    """
    model = dmax_mod.model_dmax_nm(args.pdb)
    effective = dmax_mod.resolve_dmax_nm(args.experimental, model)
    padding = dmax_mod.box_padding_nm(effective, model)
    print("MODEL_DMAX={0:.3f}".format(model))
    print("EFFECTIVE_DMAX={0:.3f}".format(effective))
    print("BOX_PADDING={0:.3f}".format(padding))
    return 0


def cmd_render_mdp(args) -> int:
    """Render an mdp template into an output file (never mutating the template).

    Replaces the ``sed -i`` in-place edit in ``run_MD.sh``. With ``--config`` the
    production overrides (nsteps/dt) are derived from the job config.
    """
    overrides = {}
    if args.config:
        overrides = mdp_mod.production_overrides(load_config(args.config))
    if args.nsteps is not None:
        overrides["nsteps"] = args.nsteps
    if args.dt is not None:
        overrides["dt"] = args.dt
    out = mdp_mod.write_rendered_mdp(args.template, args.out, overrides)
    print(out)
    return 0


def cmd_run(args) -> int:
    """Prepare the job and run (or, with --dry-run, plan) the GROMACS setup + Slurm submission.

    Always performs the local orchestration -- create the job directory tree,
    write ``job.json``, render the production mdp into the per-job mdp dir, and
    write the manifest. External steps (GROMACS setup, ``sbatch``) go through
    :class:`CommandRunner`: recorded in dry-run, executed otherwise (raising a
    clear error if the binaries are absent).
    """
    config = load_config(args.config)
    paths = _build_paths(config, args)
    runner = CommandRunner(dry_run=args.dry_run)
    status0 = manifest_mod.STATUS_PLANNED if args.dry_run else manifest_mod.STATUS_RUNNING
    manifest = manifest_mod.Manifest(config, status=status0)

    # 1. local orchestration (always real -- these are not external binaries)
    paths.create()
    os.makedirs(paths.mdp_dir, exist_ok=True)
    with open(paths.config_path, "w") as fh:
        fh.write(config.to_json())

    template_name = "lig_md.mdp" if config.system is SystemType.PROTEIN_LIGAND else "md.mdp"
    template_path = os.path.join(args.mdp_dir, template_name)
    if os.path.isfile(template_path):
        out_mdp = os.path.join(paths.mdp_dir, template_name)
        mdp_mod.write_rendered_mdp(template_path, out_mdp,
                                   mdp_mod.production_overrides(config))
        manifest.set_parameter("productionMdp", out_mdp)

    # 2. protein-ligand execution is not yet ported (planning works via 'plan')
    if config.system is SystemType.PROTEIN_LIGAND and not args.dry_run:
        manifest.set_status(manifest_mod.STATUS_FAILED).add_note(
            "protein-ligand execution not yet ported; use 'plan' or --dry-run")
        manifest.write(paths.manifest_path)
        print("protein-ligand execution is not yet ported; use --dry-run to plan.",
              file=sys.stderr)
        return 2

    # 3. box padding requires the structure
    if not (args.pdb and os.path.isfile(args.pdb)):
        manifest.set_status(manifest_mod.STATUS_FAILED).add_note("missing --pdb")
        manifest.write(paths.manifest_path)
        print("error: --pdb is required to run (needed for box padding).", file=sys.stderr)
        return 2
    model = dmax_mod.model_dmax_nm(args.pdb)
    padding = dmax_mod.box_padding_nm(dmax_mod.resolve_dmax_nm(config.dmax_nm, model), model)
    manifest.set_parameter("boxPaddingNm", padding)

    ions_mdp = os.path.join(paths.mdp_dir, "ions.mdp")
    if os.path.isfile(os.path.join(args.mdp_dir, "ions.mdp")) and not os.path.isfile(ions_mdp):
        import shutil
        shutil.copy(os.path.join(args.mdp_dir, "ions.mdp"), ions_mdp)

    try:
        # 4. GROMACS setup pipeline
        for label, argv in gmx.setup_pipeline(config, padding, ions_mdp):
            stdin = gmx.GENION_REPLACE_GROUP if label == "genion" else None
            runner.run(argv, label=label, cwd=paths.solvate_dir, stdin=stdin)
            manifest.add_step(label, argv,
                              status=(manifest_mod.STATUS_PLANNED if args.dry_run
                                      else manifest_mod.STATUS_COMPLETED))

        # 5. Slurm submission (real job-id capture when executing)
        jobs = slurm_mod.build_pipeline(config)
        symbolic = slurm_mod.symbolic_jobids(jobs)
        real_ids = {}
        for job in jobs:
            dep_ids = symbolic if args.dry_run else {
                d: real_ids.get(d, d) for d in job.depends_on}
            argv = slurm_mod.sbatch_command(
                job, args.slurm_dir, paths.simulation_dir, dep_ids, config.partition_flag())
            result = runner.run(argv, label="sbatch:" + job.name, capture=(not args.dry_run))
            if not args.dry_run and getattr(result, "stdout", None):
                real_ids[job.name] = result.stdout.decode().strip().split(";")[0]
            manifest.add_step("sbatch:" + job.name, argv,
                              status=(manifest_mod.STATUS_PLANNED if args.dry_run
                                      else manifest_mod.STATUS_COMPLETED))
    except MissingDependencyError as exc:
        manifest.set_status(manifest_mod.STATUS_FAILED).add_note(str(exc))
        manifest.write(paths.manifest_path)
        print("error: {0}".format(exc), file=sys.stderr)
        return 1

    manifest.set_status(manifest_mod.STATUS_PLANNED if args.dry_run
                        else manifest_mod.STATUS_COMPLETED)
    manifest.write(paths.manifest_path)
    print("{0}: job prepared in {1}".format(
        "planned" if args.dry_run else "submitted", paths.simulation_dir))
    print("manifest: {0}".format(paths.manifest_path))
    return 0


def cmd_analyze(args) -> int:
    """Collect SAXS fit results, write a ranked summary, and fold metrics into a manifest.

    Reads an analysis directory produced by a (real) SAXS-fitting run; needs no
    ATSAS binaries. Writes ``saxs_summary.csv`` and a manifest to the output dir.
    """
    config = load_config(args.config) if args.config else None
    out_dir = args.out_dir or args.analysis_dir
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    results = analysis_mod.collect_saxs_results(args.analysis_dir)
    csv_path = os.path.join(out_dir, "saxs_summary.csv")
    analysis_mod.write_summary_csv(results["records"], csv_path)

    manifest = analysis_mod.build_analysis_manifest(
        results, config=config, summary_csv=csv_path, saxs_dat=args.saxs)
    manifest_path = args.manifest or os.path.join(out_dir, "analysis_manifest.json")
    manifest.write(manifest_path)

    summary = results["summary"]
    print("SAXS analysis")
    print("-" * 60)
    print("  models scored : {0}".format(len(results["records"])))
    print("  best chi^2    : {0}".format(summary["bestChi2"]))
    print("  best frame    : {0}".format(summary["bestFrame"]))
    print("  mean Rg       : {0}".format(summary["rgMean"]))
    print("  summary csv   : {0}".format(csv_path))
    print("  manifest      : {0}".format(manifest_path))
    if not results["records"]:
        print("  [note] no crysol_summary.txt found in {0}".format(args.analysis_dir))
    return 0


def cmd_postprocess(args) -> int:
    """Clean each repeat's trajectory and combine (trjcat/align) for clustering.

    Builds the GROMACS command plan from the job layout and runs it through
    CommandRunner (recorded in dry-run, executed otherwise). Produces
    ``combined_aligned.xtc`` + ``final.pdb`` in the production dir, which the
    ``cluster`` stage consumes.
    """
    from . import postprocess as pp

    config = load_config(args.config)
    paths = _build_paths(config, args)
    runner = CommandRunner(dry_run=args.dry_run)
    reference_gro = os.path.join(paths.pdb2gmx_dir, "GMX.gro")
    final_xtcs = []
    try:
        for i in range(1, config.n_repeats + 1):
            proc = paths.processed_dir(i)
            if not args.dry_run:
                os.makedirs(proc, exist_ok=True)
            rep = paths.repeat_dir(i)
            src_gro = os.path.join(rep, args.deffnm + ".gro")
            src_tpr = os.path.join(rep, args.deffnm + ".tpr")
            src_xtc = os.path.join(rep, args.deffnm + ".xtc")
            for label, argv, stdin in pp.plan_repeat_processing(
                    src_gro, src_tpr, src_xtc, reference_gro, skip=args.skip):
                runner.run(argv, label="rep{0}:{1}".format(i, label), cwd=proc, stdin=stdin)
            final_xtcs.append(os.path.join(proc, pp.FINAL_XTC))

        if not args.dry_run:
            os.makedirs(paths.production_dir, exist_ok=True)
        reference_pdb = os.path.join(paths.processed_dir(1), pp.FINAL_PDB)
        for label, argv, stdin in pp.plan_combine(final_xtcs, reference_pdb, skip=args.skip):
            runner.run(argv, label=label, cwd=paths.production_dir, stdin=stdin)
        if not args.dry_run and os.path.isfile(reference_pdb):
            import shutil
            shutil.copy(reference_pdb, os.path.join(paths.production_dir, pp.FINAL_PDB))
    except MissingDependencyError as exc:
        print("error: {0}".format(exc), file=sys.stderr)
        return 1

    print("postprocess: {0} {1} command(s)".format(
        len(runner.planned), "planned" if args.dry_run else "executed"))
    print("combined trajectory: {0}".format(
        os.path.join(paths.production_dir, pp.COMBINED_ALIGNED_XTC)))
    return 0


def cmd_cluster(args) -> int:
    """Run CLoNe/PCA structural clustering on a trajectory (needs mdtraj/sklearn/scipy)."""
    from . import structural as structural_mod

    try:
        result = structural_mod.run_clustering(
            args.traj, args.topo, args.out_dir,
            at_sel=args.at_sel, pca=args.pca, pdc=args.pdc,
            n_resize=args.n_resize, filt=args.filt, verbose=args.verbose)
    except MissingDependencyError as exc:
        print("error: {0}".format(exc), file=sys.stderr)
        return 1
    print("structural clustering")
    print("-" * 60)
    print("  clusters found : {0}".format(result["nClusters"]))
    print("  summary        : {0}".format(result["outputFiles"][0]))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="automd_saxs",
        description="All-atom MD + SAXS refinement workflow (Phase 1).",
    )
    sub = parser.add_subparsers(dest="command")

    plan = sub.add_parser("plan", help="Validate config and print the planned steps (no execution).")
    plan.add_argument("--config", required=True, help="job.json or legacy configurations.txt")
    plan.add_argument("--pdb", help="input PDB, used to compute model Dmax / box padding")
    plan.add_argument("--work-dir", help="base directory for the job tree (default: cwd)")
    plan.add_argument("--slurm-dir", default="slurms", help="directory holding the .slurm scripts")
    plan.add_argument("--mdp-dir", default="mdp_files", help="directory holding the .mdp templates")
    plan.add_argument("--manifest", help="write the plan manifest JSON to this path")
    plan.add_argument("--dry-run", action="store_true", default=True,
                      help="plan only; never execute (always on for 'plan')")
    plan.set_defaults(func=cmd_plan)

    dmax = sub.add_parser("dmax", help="Print model/effective Dmax and box padding.")
    dmax.add_argument("--pdb", required=True, help="input PDB structure")
    dmax.add_argument("--experimental", help="experimental Dmax (nm) or 'Model'")
    dmax.set_defaults(func=cmd_dmax)

    render = sub.add_parser("render-mdp", help="Render an mdp template to a file (no in-place edit).")
    render.add_argument("--template", required=True, help="source .mdp template (read only)")
    render.add_argument("--out", required=True, help="output .mdp path")
    render.add_argument("--config", help="job.json / configurations.txt for production overrides")
    render.add_argument("--nsteps", type=int, help="override nsteps")
    render.add_argument("--dt", type=float, help="override dt (ps)")
    render.set_defaults(func=cmd_render_mdp)

    run = sub.add_parser("run", help="Prepare and run (or --dry-run plan) the pipeline.")
    run.add_argument("--config", required=True, help="job.json or legacy configurations.txt")
    run.add_argument("--pdb", help="input PDB (required unless --dry-run without padding)")
    run.add_argument("--work-dir", help="base directory for the job tree (default: cwd)")
    run.add_argument("--slurm-dir", default="slurms", help="directory holding the .slurm scripts")
    run.add_argument("--mdp-dir", default="mdp_files", help="directory holding the .mdp templates")
    run.add_argument("--dry-run", action="store_true", default=False,
                     help="record commands without executing GROMACS/Slurm")
    run.set_defaults(func=cmd_run)

    analyze = sub.add_parser("analyze", help="Summarise SAXS fit results into a CSV + manifest.")
    analyze.add_argument("--analysis-dir", required=True,
                         help="directory containing crysol_summary.txt / *.fit / Rg_distr*.txt")
    analyze.add_argument("--config", help="job.json / configurations.txt for manifest context")
    analyze.add_argument("--saxs", help="experimental .dat (records the ATSAS commands)")
    analyze.add_argument("--out-dir", help="where to write summary + manifest (default: analysis dir)")
    analyze.add_argument("--manifest", help="manifest output path")
    analyze.set_defaults(func=cmd_analyze)

    post = sub.add_parser("postprocess", help="Clean + combine (trjcat/align) repeat trajectories.")
    post.add_argument("--config", required=True, help="job.json or legacy configurations.txt")
    post.add_argument("--work-dir", help="base directory for the job tree (default: cwd)")
    post.add_argument("--deffnm", default="md", help="production deffnm (default: md)")
    post.add_argument("--skip", type=int, default=2, help="frame stride (legacy default: 2)")
    post.add_argument("--dry-run", action="store_true", default=False,
                      help="record commands without executing GROMACS")
    post.set_defaults(func=cmd_postprocess)

    cluster = sub.add_parser("cluster", help="CLoNe/PCA structural clustering of a trajectory.")
    cluster.add_argument("--traj", required=True, help="trajectory file (e.g. combined .xtc)")
    cluster.add_argument("--topo", required=True, help="topology file (e.g. .pdb/.gro)")
    cluster.add_argument("--out-dir", required=True, help="output directory")
    cluster.add_argument("--at-sel", default="name CA", help="atom selection (mdtraj syntax)")
    cluster.add_argument("--pca", type=int, default=0, help="number of principal components (0=off)")
    cluster.add_argument("--pdc", type=float, default=4.0, help="CLoNe density percentile parameter")
    cluster.add_argument("--n-resize", type=float, default=4.0, help="CLoNe neighbour-matrix resize")
    cluster.add_argument("--filt", type=float, default=0.1, help="CLoNe outlier filter fraction")
    cluster.add_argument("--verbose", action="store_true", help="CLoNe verbose output")
    cluster.set_defaults(func=cmd_cluster)

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        parser.print_help()
        return 2
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
