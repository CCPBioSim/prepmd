#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:39:36 2026

@author: rob
"""

import MDAnalysis as mda
from prepmd.ligand import load_universe

def get_ligand_traj(top, traj):
    u = load_universe(top)
    u2 = u.load_new(traj)
    ligand = u2.select_atoms('not protein and not water')
    ligands = ligand.split("residue")
    for ligand in ligands:
        pass #TODO