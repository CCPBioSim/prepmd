#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

# note: NOT USED, can maybe delete

from gemmi import cif

test = "/home/rob/bad_structures/protease/6xov.cif"

def get_row_item(name, row_names, row_contents):
    index = row_names.index(name)
    return row_contents[index]

def check_missing_residues(cif_file):
    doc = cif.read_file(cif_file)  # copy all the data from mmCIF file
    block = doc[0]
    a = block.find_loop("_pdbx_unobs_or_zero_occ_residues.id")
    #element = a[0]
    l = a.get_loop()
    if l is None:
        print("No residues missing.")
    if l is not None:
        print("Found "+str(l.length())+" missing residues.")
    length = l.length()
    width = l.width()
    row_names = l.tags
    to_rebuild_ids = []
    to_rebuild_res = []
    to_rebuild_asym = []
    for row in range(length):
        row_contents = []
        for col in range(width):
            row_contents.append(l[row, col])
        resname = get_row_item("_pdbx_unobs_or_zero_occ_residues.label_comp_id", row_names, row_contents)
        resid = get_row_item("_pdbx_unobs_or_zero_occ_residues.label_seq_id", row_names, row_contents)
        asym = get_row_item("_pdbx_unobs_or_zero_occ_residues.label_asym_id", row_names, row_contents)
        print("Missing residue "+resname+" at sequence ID "+resid+"")
        to_rebuild_ids.append(resid)
        to_rebuild_res.append(resname)
        to_rebuild_asym.append(asym)
    return to_rebuild_ids, to_rebuild_res, to_rebuild_asym

def get_range():
    return
        
#ids, res, asym = check_missing_residues(test)
