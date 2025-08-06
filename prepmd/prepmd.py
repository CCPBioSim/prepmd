#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import argparse
import pathlib
import os
import glob
import random
import string
import sys

from prepmd import download
from prepmd import run
from prepmd import fix
from prepmd import model

parser = argparse.ArgumentParser(prog="prepmd",
                    description="Get structures from the PDB ready for molecular dynamics runs",
                    epilog="Also, run 'prepmd test' for the test suite")
parser.add_argument("code", help="4 or 12-character PDB code")
parser.add_argument("out", help="Output filename")
parser.add_argument("-w", "--wdir", help="Working directory")
parser.add_argument("-d", "--directory",
                    help="Input directory (will automatically check for fasta sequences and structure files here)")
parser.add_argument(
    "-f", "--fasta", help="Fasta sequence to use for filling loops")
parser.add_argument("-s", "--structure",
                    help="Input structure file (pdb or mmCif)")
parser.add_argument("-a", "--alignmentout",
                    help="Alignment output file from sequences aligned for loop filling", default="alignment_out.fasta")
parser.add_argument("-fmt", "--format",
                    help="Structure format to download (only used when no structure file is provided)", default="mmCif")
parser.add_argument(
    "-q", "--quiet", help="Do not print debug info", action="store_true")
parser.add_argument("-e", "--fixstart",
                    help="Fix pdb at the end of the process, not the start", action="store_true")
parser.add_argument("-dl", "--download",
                    help="Download the sequence from UNIPROT instead of using the pdb or an external fasta file", action="store_true")
parser.add_argument("-m", "--leavemissing",
                    help="Don't restore missing atoms", action="store_true")


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def run_test(code):
    params = {"code": code,
              'outmodel': "/home/rob/temp/missingresiduestestrun/"+code+"_ready.pdb",
              "workingdir": "/home/rob/temp/missingresiduestestrun/"+code+"_"+id_generator(6)
              }
    prep(params["code"], params["outmodel"],
         params["workingdir"], structure_format="pdb")


def prep(code, outmodel, workingdir, folder=None, fastafile=None, inmodel=None,
         alignmentout="alignment_out.fasta", structure_format="mmCif",
         quiet=False, fix_after=True, download_sequence=False, fix_missing_atoms=True):

    if not os.path.isdir(workingdir):
        pathlib.Path(workingdir).mkdir(parents=True, exist_ok=True)

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
                "Mulitple structure files in folder, please specify a structure file")
        if len(models) == 0:
            raise IOError("Couldn't find input structure file in folder")

    # download structure
    if code and inmodel == None and folder == None:
        print("Downloading structure file")
        inmodel = download.get_structure(
            code, workingdir, file_format=structure_format)

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
                "Multiple sequence files found in folder, please specify a sequence file")
        fastas = None
        print("Found fasta file "+str(fastafile)+" in "+str(folder))
    if (fastafile is None and download_sequence) or (fastas is None and download_sequence):
        try:
            print("Downloading sequence file(s)")
            sequences = download.get_uniprot_sequence(code)
            # I think the problem is that the sequence downloader should only download chains that also exist in the model
            fastafile = workingdir+"/"+code+".fasta"
            with open(fastafile, "w") as file:
                file.write(sequences)
        except Exception as e:
            print(e)
            raise IOError(
                "No sequence files found in folder or specified, could not download sequence")

    if not fastafile:
        print("Fasta sequence: "+str(fastafile))

    print("Filling in "+code)
    model.fix_missing_residues(code, fastafile, alignmentout,
                               inmodel, outmodel, workingdir)
    if fix_after:
        print("Fixing PDB")
        fix.fix(outmodel, outmodel, fix_missing_atoms=fix_missing_atoms)

    print("Restoring metadata...")
    if ".pdb" in inmodel:
        fix.restore_metadata_pdb(inmodel, outmodel)
    if ".cif" in inmodel or ".mmcif" in inmodel:
        print("not implemented yet")

    print("Testing "+code)
    run.test_sim(outmodel)
    print("Wrote "+outmodel)
    
def entry_point():
    if sys.argv[1] == "test":
        tests()
    args = parser.parse_args()
    fix_after = not args.fixstart
    if args.wdir is None:
        args.wdir = args.code+"_"+id_generator(6)
    prep(args.code, args.out, os.getcwd()+os.path.sep+args.wdir,
         folder=args.directory, fastafile=args.fasta, inmodel=args.structure,
         alignmentout=args.alignmentout, structure_format=args.format,
         quiet=args.quiet, fix_after=fix_after,
         download_sequence=args.download, fix_missing_atoms=args.leavemissing)

if __name__ == "__main__":
    entry_point()

def tests():
    os.system("")
    class style():
        RED = '\033[31m'
        GREEN = '\033[32m'
        YELLOW = '\033[33m'
        BLUE = '\033[34m'
        MAGENTA = '\033[35m'
        CYAN = '\033[36m'
        WHITE = '\033[37m'
        UNDERLINE = '\033[4m'
        RESET = '\033[0m'

    testinputs = ["6xov", "9CS5", "8RM8", "8VV2", "9B8B", "8CAE", "8QZA", "8RTL", "8RTO", "9A9W"]
    results = {}
    state = 0

    for test in range(len(testinputs)):
        try:
            print(f"Testing {testinputs[test]} ({test}/{len(testinputs)})")
            run_test(testinputs[test])
            print(f"{style.GREEN}PASSED: {test} {style.RESET}")
            results[testinputs[test]] = "PASS"
        except Exception as e:
            print(f"{style.RED}FAILED: {testinputs[test]}{style.RESET}")
            results[testinputs[test]] = e
            state = 1
    print("")
    print("RESULTS:")
    for name, result in results.items():
            if result == "PASS":
                print(f"{name}: {style.GREEN}{result}{style.RESET}")
            else:
                errtype = type(result).__name__,          # TypeError
                errfile = __file__,                  # /tmp/example.py
                errline = result.__traceback__.tb_lineno  # 2
                error = errtype+" on line "+errline+" in "+errfile
                print(f"{name}: {style.RED}{error}{style.RESET}")
    for name, result in results.items():
            if result != "PASS":
                print(f"{style.YELLOW}{name} exception: {result}{style.RESET}")
                

    sys.exit(state)


# big todo: full and thorough mmcif support
