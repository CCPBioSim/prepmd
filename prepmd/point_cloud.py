#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 15:41:18 2026

@author: rob
"""

import os
from sklearn.neighbors import NearestNeighbors
import numpy as np
from prepmd.lib.icp import icp
import mrcfile
import numpy as np
import MDAnalysis as mda


def to_point_cloud(mrcdata, voxel, contour_level):
    """
    Convert an EM density map to a point cloud.
    Args:
        mrcdata - ndarray of size resolutionXresolutionXresolution
        containg MRC data.
        voxel - voxel size, from the mrcfile library, should contain three
        member variables, x, y, and z (floats), for the voxel size in those
        dimensions.
        contour_level: density above which to add a point to the point cloud,
        a float.
    Returns:
        point cloud as an ndarray with three columns (x, y, z) and a row for
        each point.
    """
    point_cloud = []
    for x in range(len(mrcdata)):
        for y in range(len(mrcdata[x])):
            for z in range(len(mrcdata[x,y])):
                if mrcdata[x,y,z] > contour_level:
                    point_cloud.append([x*voxel.x, y*voxel.y, z*voxel.z])
    return np.array(point_cloud)

def score_pdb_map(pdb, em_map, contour_level):
    """
    For a given pdb and em_map, convert them to point clouds and score their
    similarity based on the error in an alignment between two point clouds.
    Args:
        pdb - path to a pdb file, a string
        em_map - path to an an em map file for the same structure as that pdb,
        a string.
        contour_level: density above which to add a point to the point cloud,
        a float.
    Returns:
        the error, a float
    """
    with mrcfile.open(em_map) as mrc:
        vsize = mrc.voxel_size
        point_cloud = to_point_cloud(mrc.data, mrc.voxel_size, 0.01)
    u = mda.Universe(pdb)
    atoms_pointcloud = u.atoms.positions
    transformation, err, = icp(point_cloud, atoms_pointcloud)
    return np.mean(err)
        
#score = score_pdb_map('/home/rob/temp/testtemp/6XOV_PGNMLW/6XOV_fill.B99990001.pdb',
#              "/home/rob/temp/testtemp/emd_22281.map.gz",
#              0.01)
#print(score)