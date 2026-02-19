
# prepmd
[![prepmd CI](https://github.com/CCPBioSim/mdprep/actions/workflows/python-app.yml/badge.svg)](https://github.com/CCPBioSim/mdprep/actions/workflows/python-app.yml)

A utility to automatically prepare structures from the PDB for molecular dynamics simulation and perform minimisations and simple MD simulations.

## Features
* [X] Automatically download structures, sequences and metadata from the PDB and UNIPROT
* [X] Automatically fill missing loops with modeller
* [X] Automatically add missing atoms and fix non-standard residues with pdbfixer
* [X] Automatically resolve steric clashes and minimise structures
* [X] Automatically trim together structures to be the same length
* [X] Run simple MD simulations for testing, validation and minimisation
* [X] Create 'morph' trajectories with metadynamics
* [ ] Automatically propagate metadata through to finalised structure files
* [ ] AIIDA integration

## Installation
* Install [Conda](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install) (if you don't already have it)
* Clone this repo and enter the folder: `git clone https://github.com/CCPBioSim/prepmd.git && cd prepmd` 
* Run `conda env create --name prepmd --file environment.yaml && conda activate prepmd && pip install .`
* For the MODELLER part of the workflow to work, you need to get a [modeller license key](https://salilab.org/modeller/registration.html) and add it to modeller's config.py file. If you use conda, the key will be in `envs/prep/lib/modeller-10.7/modlib/modeller/config.py` relative to the path where conda is installed.
* After installing, run `pytest` to run tests.

## Preparing structures from the PDB for simulation

### Basic example: 
`prepmd 6xov 6xov_processed.pdb` will download the structure for PDB entry `6xov`, process it and write it to `6xov_processed.pdb`.
### Using a local structure file:
 `prepmd --structure 6xov_input.pdb 6xov 6xov_processed.pdb`. You still need to supply a PDB code, as the various file formats used by prepmd require one to be present.
### Generate multiple structure files:
`prepmd 6xov 6xov_processed.pdb -n 5` will generate 5 candidate structures and select the best one as determined by MODELLER's internal metrics. Alternatively, `prepmd 6xov 6xov_processed.pdb -n 5 -em 22281 --contour 0.01` will download EMD-22281, the EMDB entry associated with 6XOV, and score the generated models based on their agreement with the EM density map.
### Use refined structures from PDB-REDO:
`prepmd 1cbs 1cbs_processed.pdb --redo` will download a refined structure from PDB-REDO, if it is available. Note: not all PDB entries have corresponding PDB-REDO entries.
### Use your own alignments and sequences to fill missing loops:
By default, `prepmd` will read missing residues from the pdb/mmcif metadata, attempt to align the missing residues with the currently present residues, and then build missing loops. You can manually provide a FASTA file containing the alignment data with `--fasta`. You can also ask prepmd to get the sequence data from UNIPROT instead, with `--download`, though this is not recommended, as the raw sequence data can be different from the PDB and cause the alignment to fail.
### Other usage notes
* `prepmd` will attempt to guess the correct file format from the filenames it's given. It won't perform implicit conversions, so make sure to start and end with the same file type.
* By default, `prepmd` will leave intermediate files in a randomly-named temporary directory. You can set the name of this directory: `prepmd --wdir 6xov_temp 6xov 6xov.cif`.
* While both pdb and mmCif are supported, using the mmCif format is strongly recommended, as the pdb format has been deprecated since 2024.
* Use `prepmd --help` for a full list of parameters. 

## Running MD simulations
`runmd` can run MD simulations using OpenMM.
### A Basic Example
 `runmd structure.cif --min_out structure_minimised.cif --traj_out traj.xtc --md_steps 5000 --step 100` will minimise and run a simulation of structure.cif using OpenMM, writing a trajectory to `traj_out.xtc`, for 5000 steps, saving one trajectory frame every 100 steps.
 If you already have a minimised structure, you can skip minimisation: `runmd structure.cif --traj_out traj.xtc --md_steps 5000 --step 100 -nomin -notest`
### Explicit solvent:
`runmd structure.cif -o structure_minimised.cif --traj_out traj.xtc --md_steps 500 --step 10 -solv tip4pew` will run a simulation with the tip4pew solvent. tip3p, tip4pew and spce are supported. You can also add pressure coupling with `--pressure 1.0` (for 1 bar). By default, simulations run with an implicit solvent equivalent to AMBER's `igb=8` option.
### Force Fields:
`runmd structure.cif -o structure_minimised.cif --traj_out traj.xtc --md_steps 500 --step 50 -ff amber14` runs with amber14. charmm36, amoeba, amber14 and amber19 are available, with charmm36 being the default.
### Equilibrate side chains:
`runmd structure.cif -o structure_minimised.cif --fix_backbone -solv tip4pew --notest` will fix the backbone in place and only equilibrate side chains.
### Create a morph trajectory:
`runmd pre.cif -m post.cif -o minimised_out.pdb`  will create a trajectory that smoothly transitions between pre.cif and post.cif. This trajectory is created using OpenMM's metadynamics features. Note: this should only be used for visualisation/illustration as trajectories created this way are arbitrary representations of structural transitions that aren't guaranteed to represent the underlying physics and biology.
If you have two files for the same structure which aren't aligned (e.g. they have slightly different starting/ending residues), you can trim the ends to align them: `aligntogether pre.cif post.cif pre_cropped.cif post_cropped.cif`
### Other usage notes:
* Set the numerical integrator with the `-i` flag. This can be either `VariableLangevinIntegrator` or `LangevinMiddleIntegrator`. By default, `runmd` will attempt to use the latter, and fall back to the former if the simulation becomes numerically unstable.
* The default settings result in a rather loose coupling to the heat bath. You can change this with the `-f` or `--friction` argument, which specified the friction coefficient coupling the system to the heat bath. Running a simulation with explicit solvent will also result in tighter coupling.
* By default, `runmd` will try to select the most optimal nonbonded interaction method, but this can be overridden with `-nb` or `--nonbonded`, which can be one of `PME`, `CutoffPeriodic`, or `CutoffNonPeriodic`
* By default, `runmd` will constrain the length of all bonds involving a hydrogen atom, which can allow for longer timesteps at the cost of some accuracy. This can be disabled by setting `-c None` or `--constraints None`. This setting is also disabled if the backbone is fixed.
* Use `runmd --help` for a full list of parameters. 

### What next?
* Though you can run simple MD simulations with prepmd, for more in-depth MD we recommend using real MD software such as GROMACS, AMBER, NAMD or OpenMM.
* If you're looking to generate an atomistic structure file that matches your EM map as closely as possible, you can use a flexible fitting tool such as [TEMPy-ReFF](https://gitlab.com/topf-lab/tempy-reff).

## Licence
AGPLv3

## Contributors
prepmd is developed by Rob Welch. Thanks to Harry Swift for helping set up the CI. This project is funded by [DRIIMB](https://driimb.org/). prepmd makes use of 

## Dependencies
* OpenMM
* PDBFixer
* BioPython
* MODELLER
* pdb2pqr
* mrcfile
* icp
