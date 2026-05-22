#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integration tests
"""

import openmm.unit as unit
from prepmd.run import run
from prepmd.prep import prep
from prepmd.align_together import align_together
import os


def run_sim(code, path, fmt="pdb", ):
    prep(code,
         str(path)+os.path.sep+code+"."+fmt,
         str(path)+os.path.sep+"testout"+os.path.sep+code+"_test",
         download_format=fmt, no_modeller=True)


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
        
    def test_bestpdb(self, tmp_path):
        path = str(tmp_path)
        code = "6XOV"
        prep(code,
             str(path)+os.path.sep+code+"."+"pdb",
             str(path)+os.path.sep+"testout"+os.path.sep+code+"_test",
             download_format="pdb",
             num_models=2, em_map="22281", em_contour=0.01,
             no_modeller=True)
        
    def test_pqr(self, tmp_path):
        path = str(tmp_path)
        code = "1UBQ"
        prep(code,
             str(path)+os.path.sep+code+"."+"pdb",
             str(path)+os.path.sep+"testout"+os.path.sep+code+"_test",
             download_format="pdb",
             pqr_out=str(path)+os.path.sep+code+"."+"pqr",
             no_modeller=True)

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

    def test_minimise_run(self, tmp_path):
        run(test_file, minimised_structure_out=str(tmp_path)+"101m_prof_min.cif",
            traj_out=str(tmp_path)+"101M_proc.xtc", md_steps=100, step=25,
            solvent=None, pressure=None,
            write_params=str(tmp_path)+sep+"params.json",
            thermo_out_file=str(tmp_path)+sep+"thermo.txt",
            checkpoint_output=str(tmp_path)+sep+"checkpoint.dat")

    def test_variable_langevin(self, tmp_path):
        run(test_file,
            traj_out=str(tmp_path)+"101M_proc.xtc", md_steps=10, step=2,
            integrator_str="VariableLangevinIntegrator",
            solvent="tip4pew", pressure=1.0*unit.bar, minimise=False,
            test_run=False,
            write_params=str(tmp_path)+sep+"params.json",
            thermo_out_file=str(tmp_path)+sep+"thermo.txt",
            checkpoint_output=str(tmp_path)+sep+"checkpoint.dat")

    def test_amber14(self, tmp_path):
        run(test_file,
            traj_out=str(tmp_path)+"101M_proc.xtc", md_steps=5, step=1,
            solvent="tip3p", forcefield_str="amber14", minimise=False,
            test_run=False,
            md_timestep=0.001*unit.picoseconds,
            write_params=str(tmp_path)+sep+"params.json",
            thermo_out_file=str(tmp_path)+sep+"thermo.txt",
            checkpoint_output=str(tmp_path)+sep+"checkpoint.dat")

    def test_fix_backbone(self, tmp_path): # huge memory leak????
        run(test_file,
            traj_out=str(tmp_path)+"101M_proc.xtc", md_steps=5, step=1,
            solvent="tip4pew", pressure=1.0*unit.bar, forcefield_str="amber14",
            fix_backbone=True, minimise=True, test_run=False,
            md_timestep=0.001*unit.picoseconds,
            write_params=str(tmp_path)+sep+"params.json",
            thermo_out_file=str(tmp_path)+sep+"thermo.txt",
            checkpoint_output=str(tmp_path)+sep+"checkpoint.dat")
        
    def test_metamorph(self, tmp_path):
        testpath = os.path.dirname(file_path)+sep+"test_data"+sep
        run(testpath+"6xou_cropped.pdb",
            metadynamics_morph = testpath+"6xov_cropped.pdb",
            minimised_structure_out = str(tmp_path)+"_min.pdb",
            meta_rmsd_threshold_nm = 0.32)
        # note: the threshold being high means this isn't a very rigorous
        # test -  this test always passes on my local machine but it stalls
        # on the github CI for unknown reasons - possibly to do with the RNG?

    def test_ligand(self, tmp_path):
        testpath = os.path.dirname(file_path)+sep+"test_data"+sep
        run(testpath+"3PDT_out.pdb",
            traj_out=str(tmp_path)+"3PDT.xtc", md_steps=6, step=2,
            ligands = [testpath+"ADP.sdf"],
            write_params=str(tmp_path)+sep+"params.json",
            thermo_out_file=str(tmp_path)+sep+"thermo.txt",
            checkpoint_output=str(tmp_path)+sep+"checkpoint.dat",
            minimise=True, test_run=True)

test_aln1 = os.path.dirname(file_path)+sep+"test_data"+sep+"6xov_prep.pdb"
test_aln2 = os.path.dirname(file_path)+sep+"test_data"+sep+"6xou_prep.pdb"

# can't be tested without a modeller key because it uses modeller for alignment
#class TestAlignTogether:
#    
#    def test_align_together(self, tmp_path):
#        align_together(test_aln1, test_aln2, 
#                       str(tmp_path)+sep+"6xov_cropped.pdb",
#                       str(tmp_path)+sep+"6xou_cropped.pdb", "6xov", "6xou")