"""MD stages: minimise, equilibrate, production (OpenMM, imported lazily).

Explicit-solvent equilibrium MD (PME, LangevinMiddleIntegrator, MonteCarloBarostat
for NPT) -- the AutoMD-SAXS refinement protocol, using the OpenMM coding patterns
from the BilboMD worker (energy gating after minimisation, StateDataReporter +
DCDReporter, platform selection). OpenMM is imported only when these run.

Each stage is independently callable so the orchestration level stays testable
even though the numerical run needs OpenMM + (ideally) a GPU.
"""

import math
import os
import subprocess

from ..command_runner import MissingDependencyError
from .schema import OpenMMConfig

# Fail-fast gate for a failed minimisation. For an EXPLICIT-solvent system the
# minimised potential energy is large and NEGATIVE (the water box dominates;
# values of -1e6..-1e7 kJ/mol are normal and healthy), so we must not reject on
# magnitude. A failed minimisation instead shows up as NaN/inf or a large
# POSITIVE energy (unrelieved steric clashes). This threshold is that positive
# clash ceiling, not an absolute-value bound.
MAX_POSITIVE_ENERGY_KJ = 1.0e7


def _require_openmm():
    try:
        import openmm
        from openmm import app, unit
    except ImportError as exc:
        raise MissingDependencyError(
            "OpenMM MD requires the openmm package (conda-forge). "
            "Original import error: {0}".format(exc))
    return openmm, app, unit


def gpu_count() -> int:
    """Number of GPUs visible to this process (for parallelising repeats).

    Honours the k8s / CUDA convention: if ``CUDA_VISIBLE_DEVICES`` is set (the
    device plugin sets it per pod) its entry count wins; otherwise fall back to
    ``nvidia-smi -L``. Returns >=1 (1 means "no GPU fan-out", i.e. run repeats
    sequentially). Never raises.
    """
    cvd = os.environ.get("CUDA_VISIBLE_DEVICES")
    if cvd is not None:
        devs = [d for d in cvd.split(",") if d.strip() != ""]
        return max(1, len(devs))
    try:
        out = subprocess.run(["nvidia-smi", "-L"], capture_output=True,
                             text=True, timeout=10)
        n = len([ln for ln in out.stdout.splitlines() if ln.strip().startswith("GPU ")])
        return max(1, n)
    except Exception:  # noqa: BLE001 - GPU probing is best-effort
        return 1


def _select_platform(openmm, requested):
    if requested and requested.lower() != "auto":
        return openmm.Platform.getPlatformByName(requested)
    # fastest available: CUDA > OpenCL > CPU
    names = [openmm.Platform.getPlatform(i).getName()
             for i in range(openmm.Platform.getNumPlatforms())]
    for preferred in ("CUDA", "OpenCL", "CPU"):
        if preferred in names:
            return openmm.Platform.getPlatformByName(preferred)
    return None


def _platform_properties(config, platform, device_index):
    """Platform properties: pin to a GPU (device_index) and pick precision.

    ``device_index`` selects a specific CUDA/OpenCL device so parallel repeats
    each run on their own GPU. HMR (larger timestep) uses mixed precision for
    energy stability; otherwise the platform default (single on CUDA) is fastest.
    """
    if platform is None:
        return None
    name = platform.getName()
    props = {}
    if name == "CUDA":
        if device_index is not None:
            props["CudaDeviceIndex"] = str(device_index)
        if getattr(config, "hmr", False):
            props["CudaPrecision"] = "mixed"
    elif name == "OpenCL":
        if device_index is not None:
            props["OpenCLDeviceIndex"] = str(device_index)
        if getattr(config, "hmr", False):
            props["OpenCLPrecision"] = "mixed"
    return props or None


def _make_simulation(config, app, openmm, unit, topology, system, device_index=None):
    """Build a Simulation, pinning to ``device_index`` with the right properties."""
    platform = _select_platform(openmm, config.platform)
    props = _platform_properties(config, platform, device_index)
    if platform is not None and props:
        return app.Simulation(topology, system, _integrator(openmm, unit, config),
                              platform, props)
    return app.Simulation(topology, system, _integrator(openmm, unit, config), platform)


def build_system(config: OpenMMConfig, app, unit, topology):
    """Create an OpenMM System with PME and H-bond constraints.

    When ``config.hmr`` is set, Hydrogen Mass Repartitioning transfers mass onto
    hydrogens (``hydrogenMass=1.5 amu``) so the integration timestep can be ~2x
    larger (4 fs) for ~2x throughput at negligible accuracy cost.
    """
    forcefield = app.ForceField(*config.forcefield_files())
    kwargs = dict(
        nonbondedMethod=app.PME,
        nonbondedCutoff=config.nonbonded_cutoff_nm * unit.nanometer,
        constraints=app.HBonds,
    )
    if getattr(config, "hmr", False):
        kwargs["hydrogenMass"] = 1.5 * unit.amu
    return forcefield.createSystem(topology, **kwargs)


def _integrator(openmm, unit, config):
    integ = openmm.LangevinMiddleIntegrator(
        config.temperature_K * unit.kelvin,
        config.friction_per_ps / unit.picosecond,
        config.timestep_fs * unit.femtosecond,
    )
    if config.seed is not None:
        integ.setRandomNumberSeed(int(config.seed))
    return integ


def minimize(config: OpenMMConfig, solvated_pdb: str, out_pdb: str):
    """Energy-minimise; raise if the post-minimisation energy is unphysical."""
    openmm, app, unit = _require_openmm()
    pdb = app.PDBFile(solvated_pdb)
    system = build_system(config, app, unit, pdb.topology)
    simulation = _make_simulation(config, app, openmm, unit, pdb.topology, system)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy(maxIterations=config.minimize_max_iterations)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    # Healthy explicit-solvent minimisation yields a large NEGATIVE energy; only
    # NaN/inf or a large POSITIVE energy (unrelieved clashes) indicates failure.
    if math.isnan(energy) or math.isinf(energy) or energy > MAX_POSITIVE_ENERGY_KJ:
        raise RuntimeError("minimisation produced unphysical energy: {0} kJ/mol".format(energy))
    with open(out_pdb, "w") as handle:
        app.PDBFile.writeFile(pdb.topology, state.getPositions(), handle)
    return {"potentialEnergyKJ": energy, "output": out_pdb}


def equilibrate(config: OpenMMConfig, minimized_pdb: str, out_state: str):
    """NVT then NPT equilibration; writes a restart state (positions+velocities)."""
    openmm, app, unit = _require_openmm()
    pdb = app.PDBFile(minimized_pdb)
    system = build_system(config, app, unit, pdb.topology)
    system.addForce(openmm.MonteCarloBarostat(1.0 * unit.bar, config.temperature_K * unit.kelvin))
    simulation = _make_simulation(config, app, openmm, unit, pdb.topology, system)
    simulation.context.setPositions(pdb.positions)
    simulation.context.setVelocitiesToTemperature(config.temperature_K * unit.kelvin)
    simulation.step(config.equilibration_steps())
    simulation.saveState(out_state)
    return {"output": out_state}


def production(config: OpenMMConfig, equilibrated_state: str, minimized_pdb: str,
               out_dcd: str, log_path: str, device_index=None):
    """Run one production replicate, writing a DCD trajectory + state log.

    ``device_index`` pins this replicate to a specific GPU so repeats can run
    concurrently across multiple GPUs (see workflow fan-out).
    """
    openmm, app, unit = _require_openmm()
    pdb = app.PDBFile(minimized_pdb)
    system = build_system(config, app, unit, pdb.topology)
    system.addForce(openmm.MonteCarloBarostat(1.0 * unit.bar, config.temperature_K * unit.kelvin))
    simulation = _make_simulation(config, app, openmm, unit, pdb.topology, system,
                                  device_index=device_index)
    simulation.loadState(equilibrated_state)
    # loadState carries over the equilibration step counter; reset it so each
    # repeat's StateDataReporter (and the live per-repeat ns counter derived from
    # it) starts at step 0 / 0 ns.
    simulation.currentStep = 0
    simulation.reporters.append(app.DCDReporter(out_dcd, config.report_interval_steps))
    simulation.reporters.append(app.StateDataReporter(
        log_path, config.report_interval_steps, step=True, temperature=True,
        potentialEnergy=True, totalEnergy=True, speed=True))
    simulation.step(config.production_steps())
    return {"trajectory": out_dcd, "log": log_path}
