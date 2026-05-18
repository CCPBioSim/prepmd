#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 14:04:25 2026

@author: rob
"""

import pdb2pqr
import subprocess

def run_pdb2pqr(infile, outfile, ff="AMBER", ph=7):
    """
    Run PDB2PQR.

    Args:
        infile - path to input file, a string
        outfile - path to output PQR file, a string
        ff - force field to use for calculations and residue naming
    Returns:
        nothing, but writes 'outfile'
    """
    
    command = [infile, outfile, "--titration-state-method=propka",
               "--with-ph="+str(ph), "--ff="+ff, "--ffout="+ff, "--drop-water"]
    if pdb2pqr.__version__ < "3.7.0":
        print("WARNIONG: Old version of pdb2pqr detected, please install a "
        "version after 3.7.0. If you installed everything via conda, please  "
        "install pdb2pqr via pip (pip install pdb2pqr) as the conda package is "
        "abandoned.")
        subprocess.run( ["pdb2pqr"]+command, capture_output=True, text=True)
        return
    
    pdb2pqr.run_pdb2pqr(command)