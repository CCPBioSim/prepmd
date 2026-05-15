#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions
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
    "GLX": "Z",
}


def is_residue(resid):
    """
    Check if a string is a valid residue.
    Args:
        resid: the string
    Returns:
        a boolean that is true if the string is a valid residue
    """
    return resid in codes or list(set(resid)) == ["."]


def is_residue_sequence(sequence):
    """
    Check if a sequence only contains valid residues
    Args:
        sequence: an iterable
    Returns:
        true if every item in the iterable is a valid residue, otherwise false
    """
    for residue in sequence:
        if not is_residue(residue):
            return False
    return True


def three_to_one(resid, ignore_non_standard=False):
    """
    Convert the three-character residue representation found in SEQRES records
    to the one-character FASTA representation.
    Args:
        resid: the residue, a string
        ignore_non_standard: ignore non-standard residues instead of throwing
        an error, a bool
    Returns:
        the fasta residue, a string
    """
    if list(set(resid)) == ["."]: # hetatms
        return resid
    if is_residue(resid):
        return codes[resid]
    if ignore_non_standard:
        return ""
    else:
        raise ValueError("Residue not foudn")


def three_to_one_sequence(resids):
    """
    Convert the three-character residue representation found in SEQRES records
    to the one-character FASTA representation.
    Args:
        resid: an iterable of strings containing three-character residue codes
    Returns:
        the residue sequence in FASTA format
    """
    pdb_sequence = ""
    non_standard = []
    for resid in resids:
        try:
            pdb_sequence += three_to_one(resid)
        except ValueError:
            non_standard.append(resid)
    return pdb_sequence
