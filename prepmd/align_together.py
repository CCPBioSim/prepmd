#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trim the ends off two PDB files until their indices match.
Note: this is a big pile of spaghetti, it turns out doing this in a robust way
is actually kind of difficult. For example, I went to all the effort of getting
positions of residues relative to the start and end of each chain, rather than
measuring them from the start, only to discover that Bio.PDB throws weird
errors when you do for element in reversed(object).
"""

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser

#from Bio.PDB import StructureAlignment
import Bio.PDB
from Bio.PDB.mmcifio import MMCIFIO
import modeller
from prepmd import get_residues, model
import argparse
# from Bio import Align

"""
def align_biopython(s1, s2):

    p = PDBParser(QUIET=1)

    s1 = p.get_structure("1", s1)
    s2 = p.get_structure("2", s2)

    model_1 = s1[0]

    model_2 = s2[0]

    alignment = StructureAlignment(m1=model_1, m2=model_2)

    return alignment
"""

def align_modeller(s1, s2, code1="AAAA", code2="BBBB", outfile="alignment.out"):
    print("Getting missing residues from input file...")
    pdb1_res = get_residues.get_residues_pdb(s1, code1)
    pdb2_res = get_residues.get_residues_pdb(s2, code2)
    # pdb_fullseq = prepmd.get_residues.get_fullseq_pdb(inmodel, code)
    with open(code1+".seq", "w") as file:
        file.write("".join(pdb1_res))
    with open(code2+".seq", "w") as file:
        file.write("".join(pdb2_res))
    env = modeller.environ()
    print("Aligning sequences...")
    aln = modeller.Alignment(env)
    aln.append(file=code1+".seq", align_codes=code1)
    aln.append(file=code2+".seq", align_codes=code2)
    aln.salign()
    print("Succesfully aligned.", end="")
    aln.write(file=outfile)


# given an anlignment fasta, create an index of which sequences have gaps in which places
def get_gaps(alignment_fasta_path):
    alignment_fasta = model.fasta(alignment_fasta_path)
    sequence1 = "".join(alignment_fasta[0][1])
    sequence2 = "".join(alignment_fasta[1][1])
    sequence1 = sequence1.replace("*", "").replace("/", "")  # no asterisk
    sequence2 = sequence2.replace("*", "").replace("/", "")  # no asterisk
    gaps = ""
    for res1, res2 in zip(sequence1, sequence2):
        if res1 == "-" and res2 == "-":
            gaps += "G"
        if res1 == "-" and res2 != "-":
            gaps += "1"
        if res1 != "-" and res2 == "-":
            gaps += "2"
        else:
            gaps += "S"
    return gaps


def print_mask(mask):
    for i in mask:
        if i:
            print("O", end="")
        else:
            print("X", end="")
    print("")


def mark_missing_ends(gaps):
    dseq1 = False
    dseq2 = False
    missing_seq1 = []
    missing_seq2 = []
    for i in gaps:
        if i == "1":
            dseq1 = True
        if i == "2":
            dseq2 = True
        if i == "G":
            dseq1 = True
            dseq2 = True
        missing_seq1.append(dseq1)
        missing_seq2.append(dseq2)
    return missing_seq1, missing_seq2


# Find the innermost gap at each end at mark anything further than that gap for deletion
def mark_innermost_gap(gaps):

    gaps_half1 = gaps[:len(gaps)//2]
    gaps_half2 = gaps[len(gaps)//2:]

    delete_half1_seq1, delete_half1_seq2 = mark_missing_ends(
        ''.join(reversed(gaps_half1)))

    delete_half1_seq1.reverse()
    delete_half1_seq2.reverse()

    delete_half2_seq1, delete_half2_seq2 = mark_missing_ends(gaps_half2)

    missing_seq1 = delete_half1_seq1 + delete_half2_seq1
    missing_seq2 = delete_half1_seq2 + delete_half2_seq2

    return missing_seq1, missing_seq2

# get the exact number of residues at each end to be deleted
def get_ends_deletion(missing_seq):
    delete_start = None
    delete_end = None

    prev_res = True
    if missing_seq[0] != False:
        for res in range(len(missing_seq)):
            if missing_seq[res] != prev_res:
                delete_start = res  # indexing/ob1 error!
                break
            prev_res = missing_seq[res]
    delete_no_start = delete_start

    prev_res = True
    if missing_seq[-1] != False:
        for res in reversed(range(len(missing_seq))):
            if missing_seq[res] != prev_res:
                delete_end = res+1  # indexing/ob1 error!
                break
            prev_res = missing_seq[res]
    if type(delete_end) == int:
        delete_no_end = len(missing_seq) - delete_end
    else:
        delete_no_end = None
    return delete_no_start, delete_no_end


def remove_from_model(structure, number_start, number_end, verbose=True):
    # you need to build this index because the biopython objects don't let you
    # iterate through them backwards or do anything with the loop variables
    models = []
    for model in structure:
        chains = []
        for chain in model:
            residues = []
            for residue in chain:
                residues.append(residue.id)
            chains.append(residues)
        models.append(chains)

    number_removed = 0
    to_remove = []
    if number_end:
        for model_idx, model in zip(reversed(range(len(models))), reversed(models)):
            for chain_idx, chain in zip(reversed(range(len(model))), reversed(model)):
                for residue in reversed(chain):
                    if number_removed < number_end:
                        to_remove.append((model_idx, chain_idx, residue))
                        number_removed += 1

    if number_start:
        number_removed = 0
        for model_idx, model in zip(range(len(models)), models):
            for chain_idx, chain in zip(range(len(model)), model):
                for residue in chain:
                    if number_removed < number_start:
                        to_remove.append((model_idx, chain_idx, residue))
                        number_removed += 1

    for model, chain, res in to_remove:
        curr_chain = structure[model].child_list[chain]
        curr_chain.detach_child(res)

    print("Removed "+str(len(to_remove))+" residues.")

    return structure

def align_together(seq1path, seq2path, seq1out, seq2out,
                   code1="AAAA", code2="BBBB", verbose=True):

    align_modeller(seq1path, seq2path, outfile="alignment.out")
    gaps = get_gaps("alignment.out")
    
    if verbose:
        print("GAPS:")
        print(gaps)
        print("")
    
    missing_seq1, missing_seq2 = mark_innermost_gap(gaps)
    
    if verbose:
        print("Missing in sequence 1:")
        print_mask(missing_seq1)
        print("")
        print("Missing in sequence 2:")
        print_mask(missing_seq2)
    
    seq2_delete_start, seq2_delete_end = get_ends_deletion(missing_seq1)
    seq1_delete_start, seq1_delete_end = get_ends_deletion(missing_seq2)
    
    if verbose:
        print("Delete in sequence 1: "+str([seq1_delete_start,seq1_delete_end]))
        print("Delete in sequence 2: "+str([seq2_delete_start,seq2_delete_end]))
    
    if ".pdb" in seq1path and ".pdb" in seq2path:
        parser = PDBParser(PERMISSIVE=1)
    if ".mmcif" in seq1path+seq2path or ".mmCif" in seq1path+seq2path:
        parser = MMCIFParser(PERMISSIVE=1)
    seq1_structure = parser.get_structure(code1, seq1path)
    seq2_structure = parser.get_structure(code2, seq2path)
    
    if verbose:
        print("SEQ1 residues: "+str(len(list(seq1_structure.get_residues()))))
        print("SEQ2 residues: "+str(len(list(seq2_structure.get_residues()))))
    
    seq1_structure = remove_from_model(
        seq1_structure, seq1_delete_start, seq1_delete_end, verbose=True)
    seq2_structure = remove_from_model(
        seq2_structure, seq2_delete_start, seq2_delete_end, verbose=True)
    
    seq1_length = len(list(seq1_structure.get_residues()))
    seq2_length = len(list(seq2_structure.get_residues()))
    
    if verbose:
        print("SEQ1 residues: "+str(seq1_length))
        print("SEQ2 residues: "+str(seq2_length))
    
    if seq1_length != seq2_length:
        raise ValueError("The modified structures still don't have the same"
                         " length.")
    
    if ".mmcif" in seq1path+seq2path or ".mmCif" in seq1path+seq2path:
        io = MMCIFIO()
    if ".pdb" in seq1path and ".pdb" in seq2path:
        io = Bio.PDB.PDBIO()
    
    io.set_structure(seq1_structure)
    io.save(seq1out)
    io.set_structure(seq2_structure)
    io.save(seq2out)


def entry_point():
    "CLI entry point function. Uses sys.argv and argparse args object."

    parser = argparse.ArgumentParser(prog="aligntogether",
                                     description="Trim two structure files so "
                                     "that they are aligned and are the same "
                                     "length.")
    parser.add_argument("struct1in", help="Path to first input structure")
    parser.add_argument("struct2in", help="Path to second input structure")
    parser.add_argument("struct1out", help="Path to first output structure")
    parser.add_argument("struct2out", help="Path to second output structure")
    parser.add_argument(
        "-c1", "--code1", help="PDB code for file 1 (only used as a label)",
        default="AAAA")
    parser.add_argument(
        "-c2", "--code2", help="PDB code for file 2 (only used as a label)",
        default="BBBB")

    args = parser.parse_args()
    align_together(args.struct1in,
                   args.struct2in, 
                   args.struct1out,
                   args.struct2out,
                   args.c1,
                   args.c2)

if __name__ == "__main__":
    entry_point()