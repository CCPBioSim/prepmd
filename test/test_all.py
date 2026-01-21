#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integration tests
"""

import openmm.unit as unit
from prepmd.run import run
from prepmd.prep import prep
import os


def run_sim(code, path, fmt="pdb", ):
    prep(code,
         str(path)+os.path.sep+code+"."+fmt,
         str(path)+os.path.sep+"testout"+os.path.sep+code+"_test",
         download_format=fmt)


class TestPrep:

    def test_9CS5(self, tmp_path):
        run_sim("9CS5", tmp_path, fmt="pdb")

    def test_8CAE(self, tmp_path):
        run_sim("8CAE", tmp_path, fmt="pdb")

    def test_8QZA(self, tmp_path):
        run_sim("8QZA", tmp_path, fmt="pdb")

    def test_7IB8(self, tmp_path):
        run_sim("7IB8", tmp_path, fmt="cif")

    def test_9A9G(self, tmp_path):
        run_sim("9A9G", tmp_path, fmt="cif")

# removed: 6TY4, 6XOV, 9I3U, 8RTO (too slow!)


file_path = os.path.realpath(__file__)
sep = os.path.sep
test_file = os.path.dirname(file_path)+sep+"test_data"+sep+"101M_proc.cif"


class TestRun:

    def test_basic(self, tmp_path):
        run(test_file,
            traj_out=str(tmp_path)+"101M_proc.xtc", md_steps=100, step=25,
            solvent=None, pressure=None, minimise=False, test_run=False,
            write_params=str(tmp_path)+sep+"params.json",
            thermo_out_file=str(tmp_path)+sep+"thermo.txt",
            checkpoint_output=str(tmp_path)+sep+"checkpoint.dat")
        assert True

    def test_minimise_run(self, tmp_path):
        run(test_file, minimised_structure_out=str(tmp_path)+"101m_prof_min.cif",
            traj_out=str(tmp_path)+"101M_proc.xtc", md_steps=100, step=25,
            solvent=None, pressure=None,
            write_params=str(tmp_path)+sep+"params.json",
            thermo_out_file=str(tmp_path)+sep+"thermo.txt",
            checkpoint_output=str(tmp_path)+sep+"checkpoint.dat")
        assert True

    def test_variable_langevin(self, tmp_path):
        run(test_file,
            traj_out=str(tmp_path)+"101M_proc.xtc", md_steps=10, step=2,
            integrator="VariableLangevinIntegrator",
            solvent="tip4pew", pressure=1.0*unit.bar, minimise=False,
            test_run=False,
            write_params=str(tmp_path)+sep+"params.json",
            thermo_out_file=str(tmp_path)+sep+"thermo.txt",
            checkpoint_output=str(tmp_path)+sep+"checkpoint.dat")

    def test_amber14(self, tmp_path):
        run(test_file,
            traj_out=str(tmp_path)+"101M_proc.xtc", md_steps=5, step=1,
            solvent="tip3p", forcefield="amber14", minimise=False,
            test_run=False,
            md_timestep=0.001*unit.picoseconds,
            write_params=str(tmp_path)+sep+"params.json",
            thermo_out_file=str(tmp_path)+sep+"thermo.txt",
            checkpoint_output=str(tmp_path)+sep+"checkpoint.dat")
        assert True

    def test_fix_backbone(self, tmp_path):
        run(test_file,
            traj_out=str(tmp_path)+"101M_proc.xtc", md_steps=5, step=1,
            solvent="tip4pew", pressure=1.0*unit.bar, forcefield="amber14",
            fix_backbone=True, minimise=False, test_run=False,
            md_timestep=0.001*unit.picoseconds,
            write_params=str(tmp_path)+sep+"params.json",
            thermo_out_file=str(tmp_path)+sep+"thermo.txt",
            checkpoint_output=str(tmp_path)+sep+"checkpoint.dat")
