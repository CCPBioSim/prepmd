#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read residue information from structure files
"""
NO_MODELLER = False
try:
    from modeller import *
except:
    NO_MODELLER = True
from prepmd import util


def get_residues_pdb(pdb, code):
    """
    Get the fasta sequence of residues in the ATOM entries of a PDB or mmCif
    file.
    Args:
        pdb: path to pdb file, a string
        code: PDB code
    Returns:
        the fasta sequence as a string
    """
    if NO_MODELLER:
        raise ImportError("Can't run without MODELLER and a valid license key")
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
    """
    Get the fasta sequence of residues in the SEQRES records of a PDB/mmCif
    file.
    Args:
        pdb: path to pdb/mmcif file, a string
        code: PDB code
    Returns:
        the fasta sequence as a string
    """
    seqres = {}

    # pdb
    with open(pdb) as file:
        for line in file:
            if "SEQRES" in line:
                split = line.split()
                chain = split[2]
                if chain not in seqres.keys():
                    seqres[chain] = []
                sequence = split[4:]
                seqres[chain] += (sequence)

    # mmcif
    if seqres == {}:
        reading_seq = False
        with open(pdb) as file:
            for line in file:
                if "_entity_poly_seq" in line:
                    reading_seq = True
                if reading_seq:
                    if len(line.split()) == 4:
                        chain = line.split()[0]
                        if chain not in seqres.keys():
                            seqres[chain] = []
                        sequence = line.split()[2]
                        seqres[chain] .append(sequence)
                if line.startswith("#"):
                    reading_seq = False

    # convert to fasta
    fastas = []
    for chain, reses in seqres.items():
        fasta = ""
        for res in reses:
            if util.is_residue(res):
                fasta += util.three_to_one(res)
        fastas.append(fasta)
    fasta_joined = "/".join(fastas)
    chains = ":::::::::"
    if fastas == []:
        raise IOError("Couldn't get full sequence from contents of "+pdb+". "
                      "Does it contain a sequence?")
    return ">P1;"+code+"_fill\n"+chains+"\n"+fasta_joined+"*"
