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
from prepmd import ligand


def get_residues_pdb(pdb, code, get_hetatms=False):
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
    if get_hetatms:
        e.io.hetatm = True
    m = Model(e, file=pdb)
    aln = Alignment(e)
    aln.append_model(m, align_codes=code)
    aln.write(file=code+'.seq')
    with open(code+".seq") as file:
        original_fasta = file.readlines()
    return original_fasta


def get_fullseq_pdb(pdb, code, get_hetatms=False):
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
    hetatms_found = False
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
            if line.startswith("HET") and not hetatms_found and get_hetatms:
                hetatms_found = True

    # mmcif
    # get chain and entity order
    chains = []
    prev_chain = None
    if seqres == {}:
        with open(pdb) as file:
            for line in file:
                split = line.split()
                if "ATOM" in line and len(split) == 21:
                    chain = split[18]
                    if chain != prev_chain:
                        chains.append(chain)
                    prev_chain = chain
    
    if seqres == {}:
        reading_seq = False
        with open(pdb) as file:
            for line in file:
                if "_pdbx_poly_seq_scheme" in line:
                    reading_seq = True
                if reading_seq:
                    if len(line.split()) == 12:
                        res = line.split()[3]
                        chain = line.split()[9]
                        if chain not in seqres.keys():
                            seqres[chain] = []
                        seqres[chain].append(res)
                if line.startswith("#"):
                    reading_seq = False
                if line.startswith("HET") and not hetatms_found and get_hetatms:
                    hetatms_found = True

    """

    # mmcif
    if seqres == {}:
        reading_seq = False
        with open(pdb) as file:
            for line in file:
                if "_entity_poly_seq" in line:
                    reading_seq = True
                    mmcif = True
                if reading_seq:
                    if len(line.split()) == 4:
                        chain = line.split()[0]
                        if chain not in seqres.keys():
                            seqres[chain] = []
                        sequence = line.split()[2]
                        seqres[chain] .append(sequence)
                if line.startswith("#"):
                    reading_seq = False
                if line.startswith("HET") and not hetatms_found and get_hetatms:
                    hetatms_found = True
    
    mmcif_cross_check = False
    # mmcif cross-check
    # first, get order of chains, which for some reason is only contained in
    # the ATOM records
    # note: 5O80 shows that this is still not enough!!!
    chains = []
    prev_chain = None
    if mmcif:
        with open(pdb) as file:
            for line in file:
                split = line.split()
                if "ATOM" in line and len(split) == 21:
                    chain = split[18]
                    if chain != prev_chain:
                        chains.append(chain)
                    prev_chain = chain

    # 5O80 is five instances of the same chain and just gives a contextless sequence and five strand IDs
    # 10ED is all over the place and lists a chain ID somewhere after the sequence but it's not consistent
    # it does label them by entity ID, which is in the ATOM data - could use that?
    # no, i can't, the sequence data is labeled and even formatted differently depending on whether the file has multiple entities or different instances of the same entity
    # the only way this is all going to shake out is if I get the entity ID AND chain id of each residue,
    # then get the entity ID OR chain id of 
    
    # this is valid for 2O8A
    reading_seq = False
    fasta_lines = []
    if mmcif:
        with open(pdb) as file:
            for line in file:
                if "_entity_poly.pdbx_seq_one_letter_code" in line:
                    reading_seq = True
                if reading_seq:
                    fasta_lines.append(line.strip())
                if line.startswith("#") and fasta_lines != []:
                    reading_seq = False
                    break
        fasta_lines = "".join(fasta_lines).split(";")
        fasta_lines = list(filter(None, fasta_lines))
        fasta_lines.pop(0)
        
        
        
        
        
        seq = []
        seq_chains = []
        for line in fasta_lines:
            if "?" in line:
                seq_chains.append(line.split("?")[0])
            else:
                seq.append(line)
        seq_chains = ",".join(seq_chains).split(",")
        chain_dict = {}
        for chain, sequence in zip(seq_chains, seq):
            chain_dict[chain.strip()] = sequence
        
        fasta_joined = ""
        for chain in chains:
            #print(chain+": "+chain_dict[chain])
            try:
                fasta_joined += chain_dict[chain]
                mmcif_cross_check = True
            except KeyError:
                print("Warning: mmCif file does not provide sequence chain IDs"
                      ". If there is only one chain, this is fine, otherwise "
                      "loop building might fail.")
                break

    """

    if hetatms_found:
        # count hetatm residues and add an equivalent number of "." entries in
        # the fasta sequence
        universe = ligand.load_universe(pdb)
        ligands = universe.select_atoms('not protein and not water')
        num_hetatm_residues = len(ligands.split("residue"))
        last_key = sorted(seqres.keys())[-1]
        seqres[last_key] += [num_hetatm_residues*"."]

        
    if get_hetatms and not hetatms_found:
        raise ValueError("Was told to retrieve hetatms from "+pdb+" but none "
                         "were found.")

    # convert seqres to fasta
    fastas = []
    for chain, reses in seqres.items():
        fasta = ""
        for res in reses:
            if util.is_residue(res):
                fasta += util.three_to_one(res)
        fastas.append(fasta)
    fasta_joined = "/".join(fastas)
    

    if fastas == []:
        raise IOError("Couldn't get full sequence from contents of "+pdb+". "
                      "Does it contain a sequence?")
    chains = ":::::::::"
    return ">P1;"+code+"_fill\n"+chains+"\n"+fasta_joined+"*"
