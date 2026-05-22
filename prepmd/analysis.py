#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:39:36 2026

@author: rob
"""

import MDAnalysis as mda
from prepmd.ligand import load_universe
import numpy as np

def get_selection_traj(universe, selection):
    selection_traj = []
    for ts in universe.trajectory:
        current_pos = selection.positions
        selection_traj.append(current_pos)
    return np.array(selection_traj)


def dist(a, b):
    return ( (b[0] - a[0])**2 + (b[1] - a[1])**2 + (b[2] - a[2])**2 )**0.5


def get_ligand_centroid_traj(top, traj, ignore=["UNK"], output_trajs=True,
                             print_dists = True):
    u = load_universe(top)
    u2 = u.load_new(traj)
    ligand = u2.select_atoms('not protein and not water')
    ligands = ligand.split("residue")
    ligand_centroids = {}
    for ligand in ligands:
        res = list(set(ligand.resnames))[0]
        if res in ignore:
            continue
        all_traj = get_selection_traj(u2, ligand)
        if np.shape(all_traj)[1] > 1: # more than 1 atom
            centroid = np.mean(all_traj, axis=2)
        else:
            centroid = all_traj[0] # only 1 atom
        res_no = 2
        while True:
            if res not in ligand_centroids:
                ligand_centroids[res] = centroid
                break
            else:
                if res+"_"+str(res_no) not in ligand_centroids:
                    ligand_centroids[res+"_"+str(res_no)] = centroid
                    break
                else:
                    res_no += 1
    if output_trajs:
        for name, centroid in ligand_centroids.items():
            np.savetxt(name, centroid)
    if print_dists:
        for name, centroid in ligand_centroids.items():
            ligand_dist = dist(centroid[0], centroid[-1])
            print(name+" moved "+str(ligand_dist)+"A")
            
    return ligand_centroids



        
#top = "/home/rob/temp/testnewligand/min_3pdt.pdb"
#traj = "/home/rob/temp/testnewligand/out.xtc"
#centroids = get_ligand_centroid_traj(top, traj, output_trajs=False)
