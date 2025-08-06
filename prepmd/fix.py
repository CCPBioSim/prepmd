#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

from pdbfixer import PDBFixer
from openmm.app import PDBFile
    
def fix(pdb, out, fix_nonstandard_residues=True, fix_missing_atoms=False,
        add_missing_hydrogens=7.0, remove_heterogens=True,
        fix_missing_hydrogens=True):
    fixer = PDBFixer(filename=pdb)
    if fix_nonstandard_residues:
        print("Fixing nonstandard residues...")
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
    if remove_heterogens:
        fixer.removeHeterogens(True)
    if fix_missing_atoms:
        print("Fixing missing atoms...")
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
    if fix_missing_hydrogens:
        print("Adding missing hydrogens...")
        fixer.addMissingHydrogens(add_missing_hydrogens) #pH
    PDBFile.writeFile(fixer.topology, fixer.positions, open(out, 'w'))
    
def restore_metadata_pdb(pdb, fixed_pdb):
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
            
