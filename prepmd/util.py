#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

codes = {
    "ALA": "A",
    "ASX": "B",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "XAA": "X",
    "TYR": "Y",
    "GLX": "Z"
}


def is_residue(resid):
    return resid in codes


def is_residue_sequence(sequence):
    for residue in sequence:
        if not is_residue(residue):
            return False
    return True


def three_to_one(resid, ignore_non_standard=False):
    if is_residue(resid):
        return codes[resid]
    if ignore_non_standard:
        return ""
    else:
        raise ValueError("Residue not foudn")


def three_to_one_sequence(resids):
    pdb_sequence = ""
    non_standard = []
    for resid in resids:
        try:
            pdb_sequence += codes[resid]
        except KeyError:
            non_standard.append(resid)
    return pdb_sequence
