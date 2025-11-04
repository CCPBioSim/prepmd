# prepmd
A utility to automatically prepare structures from the PDB for molecular dynamics simulation.

## Features
* [X] Automatically download structures, sequences and metadata from the PDB and UNIPROT
* [X] Automatically fill missing loops with modeller
* [X] Automatically add missing atoms and convert non-standard residues with pdbfixer
* [X] Automatically propagate metadata through to finalised structure files
* [X] Automatically resolve steric clashes and minimise structures
* [X] A simple command-line interface for running molecular dynamics with OpenMM
* [ ] AIIDA integration

## Installation
* Install [Conda](https://conda-forge.org/download/) (if you don't already have it)
* Clone this repo and enter the folder: `git clone https://github.com/CCPBioSim/mdprep.git && cd prepmd` 
* Run `conda env create --name prepmd --file environment.yaml && conda activate prepmd && pip install .`
* For the modeller part of the workflow to work, you need to get a [modeller license key](https://salilab.org/modeller/registration.html) and add it to modeller's config.py file. If you use conda, the key will be in `envs/prepmd/lib/modeller-10.7/modlib/modeller/config.py` relative to the path where conda is installed.

## Preparing structures from the PDB for simulation
* A basic example: `prepmd 6xov 6xov_processed.pdb` will download the structure for PDB entry `6xov`, process it and write it to `6xov_processed.pdb`.
* If you already have a pdb file, you can instead run: `prepmd --structure 6xov_input.pdb 6xov 6xov_processed.pdb`. You still need to supply a PDB code, as the various file formats used by prepmd require one to be present.
* `prepmd` will attempt to guess the correct file formats from the filenames it's given. It won't perform implicit conversions, so make sure to start and end with the same file type.
* By default, `prepmd` will leave intermediate files in a randomly-named temporary directory. You can set the name of this directory: `prepmd --wdir 6xov_temp 6xov 6xov.cif`.
* By default, `prepmd` will read missing residues from the pdb/mmcif metadata, attempt to align the missing residues with the currently present residues, and then build missing loops. You can manually provide a FASTA file containing the alignment data with `--fasta`. You can also ask mdprep to get the sequence data from UNIPROT instead, with `--download`, though this is not recommended, as the raw sequence data can be very different from the PDB and cause the alignment to fail.
* Note: while both pdb and mmCif are supported, using the mmCif format is strongly recommended, as the pdb format has been deprecated since 2024.
* Use `prepmd --help` for a full list of parameters. 

## Running MD simulations
* `runmd` can run MD simulations using OpenMM.
* A basic example: `runmd structure.cif --min_out structure_minimised.cif --traj_out traj.xtc --md_steps 5000 --step 100` will minimise and run a simulation of structure.cif, writing a trajectory to `traj_out.xtc`, for 5000 steps, saving one trajectory frame every 100 steps.
* If you already have a minimised structure, you can skip minimisation: `runmd structure.cif --traj_out traj.xtc --md_steps 5000 --step 100 -nomin -notest`
* Solvate the simulation box: `runmd structure.cif -o structure_minimised.cif --traj_out traj.xtc --md_steps 500 --step 10 -solv tip4pew`. tip3p, tip4pew and spce are supported. You can also add pressure coupling with `--pressure 1.0` (for 1 bar)
* Run with different force fields: `runmd structure.cif -o structure_minimised.cif --traj_out traj.xtc --md_steps 500 --step 50 -ff amber14` runs with amber14. AMOEBA is also available, and amber19 is available if you have a recent version of OpenMM.
* Finally, you may wish to fix the backbone in place and just equilibrate the side chains: `runmd structure.cif -o structure_minimised.cif --fix_backbone -solv tip4pew --notest`
* Use `runmd --help` for a full list of parameters. 
