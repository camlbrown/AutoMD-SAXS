"""MD stages: minimise, equilibrate, production (OpenMM, imported lazily).

Explicit-solvent equilibrium MD (PME, LangevinMiddleIntegrator, MonteCarloBarostat
for NPT) -- the AutoMD-SAXS refinement protocol, using the OpenMM coding patterns
from the BilboMD worker (energy gating after minimisation, StateDataReporter +
DCDReporter, platform selection). OpenMM is imported only when these run.

Each stage is independently callable so the orchestration level stays testable
even though the numerical run needs OpenMM + (ideally) a GPU.
"""

import math

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


def build_system(config: OpenMMConfig, app, unit, topology):
    """Create an OpenMM System with PME and H-bond constraints."""
    forcefield = app.ForceField(*config.forcefield_files())
    return forcefield.createSystem(
        topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=config.nonbonded_cutoff_nm * unit.nanometer,
        constraints=app.HBonds,
    )


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
    simulation = app.Simulation(pdb.topology, system, _integrator(openmm, unit, config),
                                _select_platform(openmm, config.platform))
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
    simulation = app.Simulation(pdb.topology, system, _integrator(openmm, unit, config),
                                _select_platform(openmm, config.platform))
    simulation.context.setPositions(pdb.positions)
    simulation.context.setVelocitiesToTemperature(config.temperature_K * unit.kelvin)
    simulation.step(config.equilibration_steps())
    simulation.saveState(out_state)
    return {"output": out_state}


def production(config: OpenMMConfig, equilibrated_state: str, minimized_pdb: str,
               out_dcd: str, log_path: str):
    """Run one production replicate, writing a DCD trajectory + state log."""
    openmm, app, unit = _require_openmm()
    pdb = app.PDBFile(minimized_pdb)
    system = build_system(config, app, unit, pdb.topology)
    system.addForce(openmm.MonteCarloBarostat(1.0 * unit.bar, config.temperature_K * unit.kelvin))
    simulation = app.Simulation(pdb.topology, system, _integrator(openmm, unit, config),
                                _select_platform(openmm, config.platform))
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
