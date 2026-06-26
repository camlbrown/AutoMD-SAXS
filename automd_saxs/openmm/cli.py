"""OpenMM-branch command line: ``python -m automd_saxs.openmm ...``.

``plan`` validates the config and prints the planned stages without importing
OpenMM/FoXS -- the Phase 2 dry-run contract.
"""

import argparse
import sys

from . import workflow as workflow_mod
from .schema import OpenMMConfig


def load_config(path: str) -> OpenMMConfig:
    return OpenMMConfig.from_file(path)


def cmd_plan(args) -> int:
    config = load_config(args.config)
    stages = workflow_mod.plan_stages(config)

    out = ["AutoMD-SAXS OpenMM plan (dry run)", "=" * 60, "Configuration", "-" * 60]
    for line in config.to_json().splitlines():
        out.append("  " + line)
    out += ["", "Pipeline stages", "-" * 60]
    for i, (name, detail) in enumerate(stages, 1):
        out.append("  {0:2d}. {1:18} {2}".format(i, name, detail))

    if args.manifest:
        manifest = workflow_mod.new_manifest(config)
        for name, detail in stages:
            manifest.add_step(name)
        manifest.write(args.manifest)
        out += ["", "Manifest written to: {0}".format(args.manifest)]

    print("\n".join(out))
    return 0


def cmd_validate(args) -> int:
    try:
        load_config(args.config)
    except (ValueError, TypeError) as exc:
        print("invalid config: {0}".format(exc), file=sys.stderr)
        return 1
    print("config OK")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="automd_saxs.openmm",
        description="OpenMM/FoXS all-atom MD + SAXS refinement (Phase 2).")
    sub = parser.add_subparsers(dest="command")

    plan = sub.add_parser("plan", help="Validate config and print planned stages (no execution).")
    plan.add_argument("--config", required=True, help="job .json or .yaml")
    plan.add_argument("--manifest", help="write a planned manifest to this path")
    plan.set_defaults(func=cmd_plan)

    val = sub.add_parser("validate", help="Validate a job config and exit.")
    val.add_argument("--config", required=True, help="job .json or .yaml")
    val.set_defaults(func=cmd_validate)

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
