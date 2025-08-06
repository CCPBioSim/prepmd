#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

from prepmd import util
from modeller import *


def get_residues_pdb(pdb, code):
    log.none()
    e = Environ()
    m = Model(e, file=pdb)
    aln = Alignment(e)
    aln.append_model(m, align_codes=code)
    aln.write(file=code+'.seq')
    with open(code+".seq") as file:
        original_fasta = file.readlines()
    return original_fasta


def get_fullseq_pdb(pdb, code):
    seqres = {}
    with open(pdb) as file:
        for line in file:
            if "SEQRES" in line:
                split = line.split()
                chain = split[2]
                if chain not in seqres.keys():
                    seqres[chain] = []
                sequence = split[4:]
                seqres[chain] += (sequence)
    fastas = []
    for chain, reses in seqres.items():
        fasta = ""
        for res in reses:
            if util.is_residue(res):
                fasta += util.three_to_one(res)
        fastas.append(fasta)
    fasta_joined = "/".join(fastas)
    chains = ":::::::::"
    return ">P1;"+code+"_fill\n"+chains+"\n"+fasta_joined+"*"


def get_residues_pdb_old(pdb, with_missing=True, ignore_mutants=True):
    with open(pdb) as file:
        seqres = {}
        atomsres = {}
        prev_resno = -1
        missing_residue_remark = "nonumber"
        missing = {}
        for line in file:

            if line.startswith("SEQRES"):
                split = line.split()
                chain = split[2]
                if chain not in seqres.keys():
                    seqres[chain] = []
                sequence = split[4:]
                seqres[chain] += (sequence)

            if line.startswith("ATOM"):
                split = line.split()
                res = split[3]
                resno = split[5]
                chain = split[4]
                if len(chain) > 1:
                    chain = line[21]
                    res = split[2]
                    resno = split[4]
                    resno = ''.join(ch for ch in resno if ch.isdigit())
                if resno != prev_resno:
                    if chain not in atomsres.keys():
                        atomsres[chain] = {}
                    atomsres[chain][resno] = res
                    prev_resno = resno

            if "MISSING RESIDUES" in line:
                missing_residue_remark = line.split()[1]

            if line.startswith("REMARK "+missing_residue_remark):
                split = line.split()
                if len(split) == 5:
                    res, chain, position = split[2:5]
                    if chain not in missing.keys():
                        missing[chain] = {}
                    if ignore_mutants:
                        try:
                            int(position)
                        except:
                            continue
                    missing[chain][position] = res

        return seqres, atomsres, missing


def get_fasta_old(atomsres, missing, pdbcode, ignore_non_proteins=False):
    fasta_gaps = ""
    fasta_fill = ""
    for chain in atomsres.keys():
        if ignore_non_proteins:
            if not util.is_residue_sequence(atomsres[chain].values()):
                continue
        if chain not in missing.keys():
            # skip everything and just transcribe the sequence
            for index, res in atomsres[chain].items():
                if util.is_residue(res):
                    fasta_fill += util.three_to_one(res)
                    fasta_gaps += util.three_to_one(res)
            fasta_gaps += "/"
            fasta_fill += "/"
            continue
        # up to the highest index of both chains
        maxres = max(list(
            map(int, list(set(atomsres[chain].keys()).union(list(missing[chain].keys()))))))
        for residue_no in range(maxres):
            if str(residue_no) in missing[chain].keys():
                res = missing[chain][str(residue_no)]
                if util.is_residue(res):
                    fasta_fill += util.three_to_one(res)
                    fasta_gaps += "-"
            if str(residue_no) in atomsres[chain].keys():
                res = atomsres[chain][str(residue_no)]
                if util.is_residue(res):
                    fasta_fill += util.three_to_one(res)
                    fasta_gaps += util.three_to_one(res)
        fasta_gaps += "/"
        fasta_fill += "/"
    chains = ":::::::::"
    return ">P1;"+pdbcode+"\n"+"structureX:"+pdbcode+"::::::::"+"\n"+fasta_gaps+"*\n"+">P1;"+pdbcode+"_fill\n"+chains+"\n"+fasta_fill+"*"


def get_chains(model):
    chains = []
    with open(model) as file:
        for line in file:
            if line.startswith("ATOM"):
                chain = line.split[4]
                if chain not in chains:
                    chains.append(chain)


# sarah says to ask rene about this
# did adam hospital do something like this???
# do what chimera are doing
