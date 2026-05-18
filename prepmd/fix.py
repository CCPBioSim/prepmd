#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fix structures using PDBFixer
"""

from pdbfixer import PDBFixer
from openmm.app import PDBFile, PDBxFile


def fix_nonstandard_cif(cif):
    """
    Removed metadata written by MODELLER to cif files, which can prevent those
    files from being read by other software.
    """
    blocks = []
    curr_block = []
    with open(cif) as file:
        for line in file:
            if line.startswith("#"):
                blocks.append(curr_block)
                curr_block = []
            curr_block.append(line)
    blocks.append(curr_block)

    fixed_blocks = []
    for block in blocks:
        modeller_block = False
        for item in block:
            if "_modeller.version" in item:
                modeller_block = True
        if not modeller_block:
            fixed_blocks.append(block)

    outtext = ""
    for block in fixed_blocks:
        outtext += "".join(block)

    with open(cif, "w") as file:
        file.writelines(outtext)

# for some reason, i can't get pdbfixer to remove hetatm and unk residues?
def remove_hetatms_unk(pdb, out):
    lines = []
    with open(pdb) as file:
        for line in file:
            if line.startswith("HETATM"):
                continue
            if line.startswith("CONECT"):
                continue
            if "UNK" in line.split(" "):
                continue
            lines.append(line)
    with open(out, "w") as file:
        file.writelines(lines)


def fix(pdb, out, remove_heterogens=True, add_hydrogens=True, ph=7.0):
    """
    Fix a pdb file with pdbfixer.
    Arguments:
        pdb: path to a pdb or mmcif file, a string
        out: output file path, a string
        remove_heterogens - if true, will remove heterogens
        add_hydrogens - if true, will add hydrogens
        ph - the ph to aim for when adding hydrogens, a float
    """
    if ".cif" in pdb or ".mmcif" in pdb:
        fix_nonstandard_cif(pdb)
    fixer = PDBFixer(filename=pdb)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    if remove_heterogens:
        fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    if add_hydrogens:
        fixer.addMissingHydrogens(ph)
    if ".cif" in pdb or ".mmcif" in pdb:
        PDBxFile.writeFile(fixer.topology, fixer.positions, open(out, 'w'))
    else:
        PDBFile.writeFile(fixer.topology, fixer.positions, open(out, 'w'))


def restore_metadata_pdb(pdb, fixed_pdb):
    """
    Copy metadata from one pdb file to another. Useful as the output of
    pdbfixer and other tools often doesn't contain REMARKS and such.
    Args:
        pdb: original pdb file path containing metadata, a string
        fixed_pdb: new pdb file path, containing no metadata, a string
    Returns:
        nothing, but updates the contents of fixed_pdb
    """
    lines = []
    with open(pdb) as pdb:
        for line in pdb:
            if line.startswith("SEQRES"):
                break
            lines.append(line)
    with open(fixed_pdb) as file:
        for line in file:
            lines.append(line)
    with open(fixed_pdb, "w") as file:
        file.writelines(lines)
