#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 14:04:25 2026

@author: rob
"""

#import pdb2pqr

# steps:
    # run as normal but not including test simulation and fixing
    # apply pdb2pqr before fixing
    # then fix, and especially add missing atoms AND remove hetatms
    # then run as normal

#pdb2pqr.run_pdb2pqr("A")
# example usage: UBQ.pdb 1UBQ.pqr --titration-state-method=propka --with-ph=7 --ff=CHARMM --ffout=CHARMM