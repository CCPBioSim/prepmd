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
from pathlib import Path

from prepmd import download
from prepmd import run
from prepmd import fix
from prepmd import model
from prepmd import pdb2pqr as pqr
from prepmd import ligand

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
parser.add_argument("-q", "--quiet",
                    help="Do not print debug info", action="store_true")
parser.add_argument("-e", "--fixstart",
                    help="Fix pdb at the end of the process, not the start",
                    action="store_true")
parser.add_argument("-dl", "--download",
                    help="Download the sequence from UNIPROT instead of using the pdb or"
                    " an external fasta file", action="store_true")
parser.add_argument("-m", "--leavemissing",
                    help="Don't restore missing atoms", action="store_true")
parser.add_argument("-p", "--pqr",
                    help="PQR output filename. If set, will output a PQR "
                    "as well as a pdb.")
parser.add_argument("-ff", "--pqrfmt",
                    help="Force field to use for creating the PQR file. Can be"
                    "AMBER or LAMMPS. Will only output a PQR file if this is "
                    "set.", default="AMBER")
parser.add_argument("-r", "--redo", help="Get PDB from PDB-REDO",
                    action="store_true")
parser.add_argument("-n", "--num",
                    help="Number of models to create with MODELLER. The best "
                    "model will be selected based on the MODELLER objective "
                    "function.", type=int, default=1)
parser.add_argument("-em", "--em_map",
                    help="If multiple models are being generated, the best "
                    "model will be selected based on agreement with this "
                    "EM density map instead of the objective function.")
parser.add_argument("-c", "--contour",
                    help="Contour level for the EM map.", type=float)
parser.add_argument("-il", "--ignore_hetatms",
                    help="Delete hetatms/ligands from input file",
                    action="store_true")


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
         fix_missing_atoms=True, write_metadata="prepmeta.json", pqrff="AMBER",
         pqr_out=None, redo=False, num_models=1, em_map=None,em_contour=None,
         ignore_hetatms=False):
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
        pqrff: A string, can be AMBER or CHARMM. Force field to use for the
        PQR creation. Note: when I tested, CHARMM was buggy.
        pqr_out: output PQR filename, a string.
        redo: if True, will download PDB from PDB-REDO instead of the regular
        PDB.
        num_models: number of models to generate, an int. If >1, the best model
        will be selected based on MODELLER's internal scoring metric.
        em_map: path to an EM density map file (a string). If this is set,
        the best PDB will be picked based on similarity to the map.
        em_contour: contour level for the EM map, a float.
        ignore_hetatms: if this is set to True, hetams will be removed and not
        written to sdf files.
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

    if pqr_out and download_format == "mmCif":
        raise IOError("PQR generation requires the PDB format to be used. "
                      "Please specify an input and output PDB!")

    if not os.path.isdir(workingdir):
        pathlib.Path(workingdir).mkdir(parents=True, exist_ok=True)

    if inmodel:
        shutil.copyfile(inmodel,
                        workingdir+os.path.sep+code+"."
                        ""+download_format.replace("mmCif", "cif"))
        # note: modeller requires the filename to be the same as the pdb code
        # so here we change the filename
        inmodel = code+"."+download_format.replace("mmCif", "cif")
        
    if em_map:
        print(em_map)
        if not Path.is_file(Path(str(em_map))):
            print("Downloading EM map...")
            em_map = download.get_em_map(em_map, workingdir)
        else:
            shutil.copyfile(em_map,
                            workingdir+os.path.sep+os.path.basename(em_map))

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
            code, workingdir, file_format=download_format, redo=redo)
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

    model.fix_missing_residues(code, fastafile, alignmentout, inmodel,
                               outmodel, workingdir, num_models=num_models,
                               em_map=em_map, em_contour=em_contour,
                               keep_hetatms = not ignore_hetatms)
    
    # strip heterogens
    no_sim_output = False
    ligands = ligand.split_pdb_ligand(inmodel)
    if ligands and not ignore_hetatms:
        print("Wrote hetatms/ligands to "+", ".join(ligands))
        no_sim_output = True # don't minimise the ligandless structure without
        # the ligands!
    if ligands and ignore_hetatms:
        print("Ligands found but will be ignored.")
    if not ligands:
        print("No ligands found.")
    
    if fix_after:
        print("Fixing PDB")
        fix.fix(outmodel, outmodel, fix_missing_atoms=fix_missing_atoms)

    # TODO: why does this output broken pdbs?
    # also - this was due for a refactor anyway
    # maybe dump the data into another file instead
    #print("Restoring metadata...")
    #if ".pdb" in inmodel:
    #    fix.restore_metadata_pdb(inmodel, outmodel)
    #if ".cif" in inmodel or ".mmcif" in inmodel:
    #    print("Metadata restoration not implemented for mmCif (yet)")
    
    if pqr_out:
        print("Creating PQR...")
        pqr.run_pdb2pqr(outmodel, pqr_out, ff=pqrff)        

    print("Simulating "+code)
    if pqr_out:
        #run.test_sim(pqr_out)
        run.run(pqr_out, minimised_structure_out=None, md_steps=None,
            integrator_str="LangevinMiddleIntegrator",
            pressure=None, test_sim_steps=250)
    else:
        run.test_sim(outmodel, no_output = no_sim_output)
    print("Done.")
    if ligand and no_sim_output:
        print("Note: "+outmodel+" has not been minimsed. This is because there"
              " are possible ligands in the input file, and you may want to "
              "add these back to the structure before running a simulation.")
    
    with open(write_metadata, "w") as file:
        file.write(locals_json)

    if not os.path.isabs(outmodel) and not pqr_out:
        shutil.copyfile(outmodel, run_dir+os.path.sep+outmodel)
        
    if not os.path.isabs(outmodel) and pqr_out:
        shutil.copyfile(pqr_out, run_dir+os.path.sep+pqr_out)
        print("Wrote final PQR file to "+pqr_out+".")
    
    if not os.path.isabs(outmodel):
        shutil.copyfile(write_metadata, run_dir+os.path.sep+write_metadata)
    
    if not os.path.isabs(outmodel) and ligands and not ignore_hetatms:
        for lig in ligands:
            shutil.copyfile(lig, run_dir+os.path.sep+lig)


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
         download_sequence=args.download, fix_missing_atoms=args.leavemissing,
         pqrff=args.pqrfmt, pqr_out=args.pqr,
         redo=args.redo, num_models=args.num, em_map=args.em_map,
         em_contour=args.contour, ignore_hetatms=args.ignore_hetatms)


if __name__ == "__main__":
    entry_point()
