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

# will be removed in the next build
# it turns out a lot of these features are interdependent so if you set some
# flags and not others, pdbfixer will just fail in weird ways
"""
def fix2(pdb, out, fix_nonstandard_residues=True, fix_missing_atoms=False,
        add_missing_hydrogens=7.0, remove_heterogens=True,
        fix_missing_hydrogens=True, fix_missing_residues=True):
    
    Use PDBFixer to fix a PDB (or mmCif) structure file.
    Args:
        pdb: input file path, a string
        out: output file path, a string
        fix_nonstandard_residues: whether to fix nonstandard residues, a bool
        fix_misisng_atoms: whether to restore missing atoms, a bool
        add_missing_hydrogens: desired pH, a float
        remove_heterogens: whether to remove heterogens from the file, a bool
        fix_missing_hydrogens: whether to add missing hydrogens, a bool
        
    if fix_missing_residues:
        fix_missing_atoms = True
    if ".cif" in pdb or ".mmcif" in pdb:
        fix_nonstandard_cif(pdb)
    fixer = PDBFixer(filename=pdb)
    if fix_nonstandard_residues:
        print("Fixing nonstandard residues...")
        if fix_missing_residues:
            print("Identifying missing residues...")
            fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
    if remove_heterogens:
        print("Removing heterogens...")
        fixer.removeHeterogens(False)
    if fix_missing_atoms:
        fixer.findMissingAtoms()
        print("Adding missing atoms...")
        fixer.addMissingAtoms()
    if fix_missing_hydrogens:
        print("Adding missing hydrogens...")
        fixer.addMissingHydrogens(add_missing_hydrogens)  # pH
    if ".cif" in pdb or ".mmcif" in pdb:
        PDBxFile.writeFile(fixer.topology, fixer.positions, open(out, 'w'))
    else:
        PDBFile.writeFile(fixer.topology, fixer.positions, open(out, 'w'))
"""

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
