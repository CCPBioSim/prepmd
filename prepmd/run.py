#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import time
import json
import sys


def get_prop():
    numgpus = 1
    cuda_device_index = ",".join(map(str, range((numgpus))))
    platform = Platform.getPlatformByName('CUDA')
    prop = dict(CudaPrecision='mixed', CudaDeviceIndex=cuda_device_index) if platform.getName(
    ) == 'CUDA' else dict()
    return platform, prop


def variable_minimise(pdb, out, max_iterations=1000, timestep=0.002*picoseconds):
    # platform, prop = get_prop()
    pdb = PDBFile(pdb)
    forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
                                     nonbondedCutoff=1*nanometer, constraints=HBonds)
    integrator = VariableLangevinIntegrator(300*kelvin, 1/picosecond, 0.001)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.reporters.append(PDBReporter('output.pdb', 5000))
    simulation.minimizeEnergy(maxIterations=max_iterations)
    PDBFile.writeModel(pdb.topology, simulation.context.getState(
        getPositions=True).getPositions(asNumpy=True), file=open(out, "w"), keepIds=True)


def test_sim(pdb, minimise=10, steps=10, timestep=0.002*picoseconds):
    # platform, prop = get_prop()
    pdb = PDBFile(pdb)
    forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
                                     nonbondedCutoff=1*nanometer, constraints=HBonds)
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, timestep)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy(maxIterations=minimise)
    simulation.step(steps)
