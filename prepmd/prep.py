#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prepare structures from the PDB for molecular dynamics.
"""

import argparse
import pathlib
import os
import glob
import random
import string
import shutil
import copy
import json

from prepmd import download
from prepmd import run
from prepmd import fix
from prepmd import model

parser = argparse.ArgumentParser(prog="prepmd",
                                 description="Get structures from the PDB ready for "
                                 "molecular dynamics runs")
parser.add_argument("code", help="4 or 12-character PDB code")
parser.add_argument("out", help="Output filename")
parser.add_argument("-w", "--wdir", help="Working directory")
parser.add_argument("-d", "--directory",
                    help="Input directory (will automatically check for fasta"
                    " sequences and structure files here)")
parser.add_argument(
    "-f", "--fasta", help="Fasta sequence to use for filling loops")
parser.add_argument("-s", "--structure",
                    help="Input structure file (pdb or mmCif)")
parser.add_argument("-a", "--alignmentout",
                    help="Alignment output file from sequences aligned for loop filling",
                    default="alignment_out.fasta")
parser.add_argument("-fmt", "--dlformat",
                    help="Structure format to download (only used when no structure file is "
                    "provided)",
                    default=None)
parser.add_argument(
    "-q", "--quiet", help="Do not print debug info", action="store_true")
parser.add_argument("-e", "--fixstart",
                    help="Fix pdb at the end of the process, not the start",
                    action="store_true")
parser.add_argument("-dl", "--download",
                    help="Download the sequence from UNIPROT instead of using the pdb or"
                    " an external fasta file", action="store_true")
parser.add_argument("-m", "--leavemissing",
                    help="Don't restore missing atoms", action="store_true")


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    """
    Generate a random 6-character ID. Used for naming scratch directories.
    Args:
        size: length of string, an int
        chars; chars to use, a string
    Returns:
        the string
    """
    return ''.join(random.choice(chars) for _ in range(size))


def prep(code, outmodel, workingdir, folder=None, fastafile=None, inmodel=None,
         alignmentout="alignment_out.fasta", download_format="mmCif",
         quiet=False, fix_after=True, download_sequence=False,
         fix_missing_atoms=True, write_metadata="prepmeta.json"):
    """
    Prepare a PDB/MMCIF structure file for simulation.
    Args:
        code: pdb code, a string
        outmodel: path to output file, a string
        workingdir: path to working directory
        folder: input folder, containing structure file/fasta sequences. Used
        instead of the pdb code if provided. A string
        fastafile: path to fasta-formatted text file containing the original
        sequence for the structure file, a string. If no fastafile is provided
        the sequence will instead be found from the structure file metadata.
        inmodel: path to input structure file, a string.
        alignmentout: output file from the alignment of the sequence from the
        structure file and the original sequence.
        download_format: format of the structure file, a string. Can be
        'mmCif' or 'pdb'.
        quiet: if true, don't print anything. boolean
        fix_after: if true, fix the pdb file after loop building. If false,
        fix it before. a boolean
        download_sequence: get the fasta sequence for the structure from
        UNIPROT instead of the pdb metadata or an external file. Note: the
        UNIPROT sequence is normally very different from the pdb sequence, so
        the alignment might fail. a boolean
        fix_missing_atoms: whether to add missing atoms with pdbfixer. A bool
    Returns:
        nothing, but writes out a file to outmodel.
    """

    # don't look at this
    locals_copy = copy.copy(locals())
    locals_json = json.dumps(locals_copy)
    
    # infer download format from output format
    if not download_format and not inmodel:
        if (".pdb") in outmodel:
            download_format = "pdb"
            print("No download format specified, downloading PDB.")
        elif (".cif") in outmodel or ".mmcif" in outmodel or (
                ".pdbx") in outmodel:
            download_format = "mmCif"
            print("No download format specified, downloading mmCif.")
    
    #infer download format from input format
    if not download_format and inmodel:
        if (".pdb") in inmodel:
            download_format = "pdb"
        if ".cif" in inmodel or ".pdbx" in inmodel or ".mmcif" in inmodel:
            download_format = "mmCif"

        

    # check that input/output are specified in the same format
    # i'm not against converting the files but it shouldn't happen implicitly
    def in_string(substr, text): return text == None or substr in text.lower()
    if in_string(".pdb", inmodel) and in_string(
            ".pdb", outmodel) and in_string("pdb", download_format):
        pass
    elif in_string(".cif", inmodel) and in_string(
            ".cif", outmodel) and in_string("cif", download_format):
        pass
    else:
        raise IOError("Inputs and outputs are in different formats! Please "
                      "use only one format (ideally cif)")

    if not os.path.isdir(workingdir):
        pathlib.Path(workingdir).mkdir(parents=True, exist_ok=True)

    if inmodel:
        shutil.copyfile(inmodel, workingdir+os.path.sep+code+"."+download_format)
        # note: modeller requires the filename to be the same as the pdb code
        # so here we change the filename
        inmodel = code+"."+download_format.replace("mmCif", "cif")

    run_dir = os.getcwd()
    os.chdir(workingdir)

    # check folder for strucutre/sequence files
    if inmodel is None and folder:
        pdbs = glob.glob(folder+'/*.pdb')
        cifs = glob.glob(folder+'/*.cif')
        models = pdbs + cifs
        if len(cifs) == 1:
            inmodel = cifs[0]
        if len(pdbs) == 1:
            inmodel = pdbs[0]
        print("Found input model:"+str(inmodel))
        if len(pdbs) > 1 or len(cifs) > 1:
            raise IOError(
                "Mulitple structure files in folder, please specify a"
                " structure file")
        if len(models) == 0:
            raise IOError("Couldn't find input structure file in folder")

    # download structure
    if code and inmodel == None and folder == None:
        print("Downloading structure file")
        inmodel = download.get_structure(
            code, workingdir, file_format=download_format)
        print("Downloaded to "+inmodel)

    # fix
    if not fix_after:
        print("Fixing structure file...")
        fix.fix(inmodel, inmodel, fix_missing_atoms=fix_missing_atoms)

    # get fasta sequence
    fastas = None
    if fastafile is None and folder:
        fastas = glob.glob(folder+'/*.fasta')
        if len(fastas) == 1:
            fastafile = fastas[0]
        if len(fastas) > 1:
            raise IOError(
                "Multiple sequence files found in folder, please specify a"
                " sequence file")
        fastas = None
        print("Found fasta file "+str(fastafile)+" in "+str(folder))
    if (fastafile is None and download_sequence) or (
            fastas is None and download_sequence):
        try:
            print("Downloading sequence file(s)")
            sequences = download.get_uniprot_sequence(code)
            fastafile = workingdir+"/"+code+".fasta"
            with open(fastafile, "w") as file:
                file.write(sequences)
        except Exception as e:
            print(e)
            raise IOError(
                "No sequence files found in folder or specified, could not"
                "download sequence")

    model.fix_missing_residues(code, fastafile, alignmentout,
                               inmodel, outmodel, workingdir)
    if fix_after:
        print("Fixing PDB")
        fix.fix(outmodel, outmodel, fix_missing_atoms=fix_missing_atoms)

    print("Restoring metadata...")
    if ".pdb" in inmodel:
        fix.restore_metadata_pdb(inmodel, outmodel)
    if ".cif" in inmodel or ".mmcif" in inmodel:
        print("Metadata restoration not implemented for mmCif (yet)")

    print("Simulating "+code)
    run.test_sim(outmodel)
    print("Done.")
    
    with open(write_metadata, "w") as file:
        file.write(locals_json)

    if not os.path.isabs(outmodel):
        shutil.copyfile(outmodel, run_dir+os.path.sep+outmodel)
        shutil.copyfile(write_metadata, run_dir+os.path.sep+write_metadata)


def entry_point():
    "CLI entry point function. Uses sys.argv and argparse args object."
    args = parser.parse_args()
    fix_after = not args.fixstart
    if args.wdir is None:
        args.wdir = args.code+"_"+id_generator(6)
    prep(args.code, args.out, os.getcwd()+os.path.sep+args.wdir,
         folder=args.directory, fastafile=args.fasta, inmodel=args.structure,
         alignmentout=args.alignmentout, download_format=args.dlformat,
         quiet=args.quiet, fix_after=fix_after,
         download_sequence=args.download, fix_missing_atoms=args.leavemissing)


if __name__ == "__main__":
    entry_point()
