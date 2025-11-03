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
#import sys
import shutil

from prepmd import download
from prepmd import run
from prepmd import fix
from prepmd import model

parser = argparse.ArgumentParser(prog="prepmd",
                    description="Get structures from the PDB ready for "
                    "molecular dynamics runs")
#                    epilog="Also, run 'prepmd test' for the test suite")
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
         fix_missing_atoms=True):
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

    # infer download format from output format
    if not download_format:
        if (".pdb") in outmodel:
            download_format = "pdb"
            print("No download format specified, downloading PDB.")
        elif (".cif") in outmodel or ".mmcif" in outmodel or (
                ".pdbx") in outmodel:
            download_format = "mmCif"
            print("No download format specified, downloading mmCif.")
    
    # check that input/output are specified in the same format
    # i'm not against converting the files but it shouldn't happen implicitly
    in_string = lambda substr, text : text == None or substr in text.lower()
    if in_string(".pdb", inmodel) and in_string(
            ".pdb", outmodel) and in_string("pdb", download_format):
        pass
    elif in_string(".cif", inmodel) and in_string(
            ".cif", outmodel) and in_string("cif", download_format):
        pass
    else:
        raise IOError("Inputs and outputs are in different formats! Please "
                      "use only one format (ideally cif)" )

    if not os.path.isdir(workingdir):
        pathlib.Path(workingdir).mkdir(parents=True, exist_ok=True)

    if inmodel:
        # todo: copy input file to working directory
        shutil.copyfile(inmodel, workingdir+os.path.sep+inmodel)

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
    
    if not os.path.isabs(outmodel):
        shutil.copyfile(outmodel, run_dir+os.path.sep+outmodel)
    
def entry_point():
    "CLI entry point function. Uses sys.argv and argparse args object."
#    if len(sys.argv) == 2:
#        if sys.argv[1] == "test":
#            tests()
    args = parser.parse_args()
    fix_after = not args.fixstart
    if args.wdir is None:
        args.wdir = args.code+"_"+id_generator(6)
    prep(args.code, args.out, os.getcwd()+os.path.sep+args.wdir,
         folder=args.directory, fastafile=args.fasta, inmodel=args.structure,
         alignmentout=args.alignmentout, download_format=args.dlformat,
         quiet=args.quiet, fix_after=fix_after,
         download_sequence=args.download, fix_missing_atoms=args.leavemissing)


def test_mmcif_support():
    genid = id_generator(6)
    prep("6xov",
         os.getcwd()+os.path.sep+"6xov"+"_"+genid+".cif",
         os.getcwd()+os.path.sep+"testout"+os.path.sep+"6xov"+"_"+genid,
         download_format="mmCif")

def test_minimise():
    genid = id_generator(6)
    prep("6TY4",
         os.getcwd()+os.path.sep+"6TY4"+"_"+genid+".cif",
         os.getcwd()+os.path.sep+"testout"+os.path.sep+"6TY4"+"_"+genid)


# note: deprecated now that ctest support has been added
# will be removed!!!

# def tests():
#     os.system("")
#     class style():
#         RED = '\033[31m'
#         GREEN = '\033[32m'
#         YELLOW = '\033[33m'
#         BLUE = '\033[34m'
#         MAGENTA = '\033[35m'
#         CYAN = '\033[36m'
#         WHITE = '\033[37m'
#         UNDERLINE = '\033[4m'
#         RESET = '\033[0m'
    
#     tests = [
#        # {"id": "6TY4", "format":"pdb"},
#        # {"id": "6XOV", "format":"pdb"},
#         {"id": "9CS5", "format":"pdb"},
#         {"id": "8CAE", "format":"pdb"},
#         {"id": "8QZA", "format":"pdb"},
#         {"id": "8RTO", "format":"mmCif"},
#         {"id": "7IB8", "format":"mmCif"},
#         {"id": "9A9G", "format":"mmCif"},
#         {"id": "9I3U", "format":"pdb"},
#         #test_minimise,
#         #test_mmcif_support
#     ]
    
#     types = {"mmCif":"cif", "cif":"cif", "pdb":"pdb"}
        
#     results = {}
#     state = 0
#     cwd = os.getcwd()

#     for test in range(len(tests)):
#         try:
#             os.chdir(cwd)
#             code = tests[test]["id"]
#             curr_format =  tests[test]["format"]
#             print(f"Testing {code} ({test}/{len(tests)})")
#             genid = id_generator(6)
#             pathlib.Path(os.getcwd()+os.path.sep+"testout").mkdir(
#                 parents=True, exist_ok=True)
#             if type(tests[test]) == dict:
#                 prep(code,
#                      os.getcwd()+os.path.sep+code+"_"+genid+"."+types[curr_format],
#                      os.getcwd()+os.path.sep+"testout"+os.path.sep+code+"_"+genid,
#                      download_format=curr_format)
#             elif callable(tests[test]):
#                 test()
#             print(f"{style.GREEN}PASSED: {test} {style.RESET}")
#             results[code] = "PASS"
#         except Exception as e:
#             print(f"{style.RED}FAILED: {test}{style.RESET}")
#             results[code] = e
#             state = 1
#     print("")
#     print("RESULTS:")
#     for name, result in results.items():
#             if result == "PASS":
#                 print(f"{name}: {style.GREEN}{result}{style.RESET}")
#             else:
#                 errtype = type(result).__name__,          # TypeError
#                 errfile = __file__,                  # /tmp/example.py
#                 errline = result.__traceback__.tb_lineno  # 2
#                 error = str(errtype)+" on line " + \
#                     str(errline)+" in "+str(errfile)
#                 print(f"{name}: {style.RED}{error}{style.RESET}")
#     for name, result in results.items():
#             if result != "PASS":
#                 print(f"{style.YELLOW}{name} exception: {result}{style.RESET}")
                
#     sys.exit(state)

if __name__ == "__main__":
    entry_point()

# big todo: full and thorough mmcif support
