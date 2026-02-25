#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 14:04:25 2026

@author: rob
"""

import pdb2pqr

def run_pdb2pqr(infile, outfile, ff="AMBER"):
    """
    Run PDB2PQR.

    Args:
        infile - path to input file, a string
        outfile - path to output PQR file, a string
        ff - force field to use for calculations and residue naming
    Returns:
        nothing, but writes 'outfile'
    """
    pdb2pqr.run_pdb2pqr([infile, outfile, "--titration-state-method=propka",
                        "--with-ph=7", "--ff="+ff, "--ffout="+ff, "--drop-water"])