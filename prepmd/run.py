#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Use openmm to run simulations
"""

from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import time
import json
import sys


def get_prop():
    """
    Get platform info used by OpenMM.
    Args:
        none.
    Returns:
        platform, prop used when initialising openMM simulations
    """
    numgpus = 1
    cuda_device_index = ",".join(map(str, range((numgpus))))
    platform = Platform.getPlatformByName('CUDA')
    prop = dict(CudaPrecision='mixed',
                CudaDeviceIndex=cuda_device_index) if platform.getName(
    ) == 'CUDA' else dict()
    return platform, prop


def variable_minimise(pdb, out, max_iterations=1000, errortol=0.001):
    """
    Minimise structure with openmm, with a variable langevin integrator.
    Args:
        pdb: input structure file (pdb or mmcif), a string
        out: path to output file, a string
        max_iterations: max iterations of the miminisation to perform, an int
        errortol: arbitrary openmm value
    """
    # platform, prop = get_prop()
    if ".cif" in pdb or ".mmcif" in pdb:
        pdb = PDBxFile(pdb)
    else:
        pdb = PDBFile(pdb)
    forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')
    system = forcefield.createSystem(pdb.topology,
                                     nonbondedMethod=app.NoCutoff,
                                     nonbondedCutoff=1*nanometer,
                                     constraints=HBonds)
    integrator = VariableLangevinIntegrator(300*kelvin, 1/picosecond, errortol)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.reporters.append(PDBReporter('output.pdb', 5000))
    simulation.minimizeEnergy(maxIterations=max_iterations)
    PDBFile.writeModel(pdb.topology, simulation.context.getState(
        getPositions=True).getPositions(asNumpy=True),
        file=open(out, "w"), keepIds=True)


def test_sim(pdb, minimise=10, steps=10, timestep=0.002*picoseconds):
    """
    Run a short test simulation with openmm.
    Args:
        pdb: input structure file (pdb or mmcif), a string
        minimise: number of minimisation seps, an int
        steps: number of MD steps to run, an int
        timestep: timestep in openmm units
    Returns:
        nothing. The only thing it might do is throw an error!
    """
    # platform, prop = get_prop()
    if ".cif" in pdb or ".mmcif" in pdb:
        pdb = PDBxFile(pdb)
    else:
        pdb = PDBFile(pdb)
    forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')
    system = forcefield.createSystem(pdb.topology,
                                     nonbondedMethod=app.NoCutoff,
                                     nonbondedCutoff=1*nanometer,
                                     constraints=HBonds)
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, timestep)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy(maxIterations=minimise)
    simulation.step(steps)
