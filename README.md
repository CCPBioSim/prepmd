
# prepmd
[![prepmd CI](https://github.com/CCPBioSim/mdprep/actions/workflows/python-app.yml/badge.svg)](https://github.com/CCPBioSim/mdprep/actions/workflows/python-app.yml)

A utility to automatically prepare structures from the PDB for molecular dynamics simulation and perform minimisations and simple MD simulations.

## Features
* Automatically download structures, sequences and metadata from the PDB, PDB-REDO, EMDB and UNIPROT
* Automatically fill missing loops with MODELLER
* Automatically add missing atoms and fix non-standard residues with pdbfixer
* Automatically resolve steric clashes and minimise structures
* Automatically align and trim together structures to be the same length
* Automatically extract and prepare hetatms\ligands for simulation
* Easily run simple MD simulations for testing, validation and minimisation
* Create 'morph' trajectories with metadynamics
* Coming soon: integration with other MD\EM workflows!

## Installation

* Install [Conda](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install) (if you don't already have it). If you have an existing conda install, make sure it can install packages from conda-forge.
* Recommended: create a new virtual environment: `conda env create --name prepmd && conda activate prepmd`
* Install prepmd from the CCPBioSim conda channel: `conda install -c CCPBioSim prepmd`
* Add your [modeller license key](https://salilab.org/modeller/registration.html) by running `prep-license <your license key>` 

## Quickstart

`prepmd 6xov 6xov_processed.cif` will download the structure for PDB entry `6xov`, process it and write it to `6xov_processed.cif`.

`runmd 6xov_processed.cif --traj_out traj.xtc --md_steps 5000` will minimise and run a simulation of 6xov_processed.cif writing a trajectory to `traj_out.xtc`, for 5000 steps. By default, `runmd` uses a minimal set of simulation parameters, which aren't likely to be accurate - check the `runmd` section of this documentation for more options

## Preparing structure files for simulation with prepmd

Note: When running `prepmd`, we recommend using .mmCif files where possible. The .pdb file format is deprecated, and is provided for legacy compatibility.

### prepmd workflow
Steps in the `prepmd` workflow:
* The structure file(s) are downloaded (if not supplied) into a working directory. `prepmd` automatically infers the file format from the file extension of the input/output files.
* `prepmd` extracts the sequence from the residues in the PDB directly and compares them to a reference sequence. By default this is the sequence described in the SEQRES entries of the structure file. The two sequences are alligned and [MODELLER](https://salilab.org/modeller/) is used to fill in the missing residues.
* Optionally, multiple models can be created, and scored based on MODELLER's internal metrics or their similarity to a reference EM density map.
* Depending on settings, HETATMS may be extracted from the structure file and saved to .sdf files. [rdkit](https://www.rdkit.org/) is used to add hydrogens and correct the geometry of the ligands.
* [PDBFixer](https://github.com/openmm/pdbfixer) is used to add missing hydrogens and remove nonstandard residues.
* Optionally, at this point, a PQR file can be output using [PDB2PQR](https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/apbs/pdb2pqr.html).
* Finally, [OpenMM](https://openmm.org/) is used to perform a test minimisation and simulation. This step ensures that the resulting file is ready for simulation and that there are no steric clashes. If the minimisation or test simulation fails, it will be retried with OpenMM's variable langevin integrator. In testing, this has successfully minimised structure files with high clash scores.
* The final, mimimised structure file will be written out. Note: if ligands are present, the non-minimised structure will be written instead - this is to allow the user to choose which ligand files to include in their final structure, which can be minimised using `runmd`.

### prepmd command-line reference
* Use `prepmd --help` for a full list of parameters. 

### Worked examples

#### Using a local structure file
 `prepmd --structure 6xov_input.pdb 6xov 6xov_processed.pdb`. You still need to supply a PDB code, as some of the file formats used by prepmd require one to be present. The code doesn't have to be a 'real' PDB code, e.g. 'AAAA' will work fine. When using this setting, the input and output files must be in the same format - prepmd doesn't perform implicit conversions!
#### Generate multiple structure files
`prepmd 6xov 6xov_processed.pdb -n 5` will generate 5 candidate structures and select the best one as determined by MODELLER's internal metrics. Alternatively, `prepmd 6xov 6xov_processed.pdb -n 5 -em 22281 --contour 0.01` will download EMD-22281, the EMDB entry associated with 6XOV, and score the generated models based on their agreement with the EM density map (using the iterative closest point algorithm). The -em setting can also point to a map file.
#### Use refined structures from PDB-REDO
`prepmd 1cbs 1cbs_processed.pdb --redo` will download a refined structure from PDB-REDO, if it is available. Note: not all PDB entries have corresponding PDB-REDO entries.
#### Use your own alignments and sequences to fill missing loops
By default, `prepmd` will read missing residues from the pdb/mmcif SEQRES records, attempt to align the missing residues with the currently present residues, and then build missing loops with MODELLER. You can manually provide an aligned FASTA file containing the the complete and incomplete sequences with `--fasta`.  You can also ask prepmd to get the sequence data from UNIPROT instead, with `--download`, though this is not recommended, as the raw sequence data can be substantially different from the PDB and cause the alignment to fail.
#### Handling ligands
By default, `prepmd` strips all HETATMS records from structure files. To keep these files, run prepmd with `--hetatms`, which will strip the HETATMS records from the structure file and write each residue to its own SDF file. The co-ordinates inside the SDF files correspond to the co-ordinates of the ligands in the structure file, so the ligands can be added back into the original structure easily. `prepmd` uses [rdkit](https://www.rdkit.org/) to add hydrogens and correct the geometry of small molecules. 
#### Working directory
By default, `prepmd` will leave intermediate files in a randomly-named temporary directory. You can set the name of this directory: `prepmd --wdir 6xov_temp 6xov 6xov.cif`.
#### Other notes
Warning: `prepmd`'s output does not contain any of the metadata found in the original pdb. This is an intentional omission - a lot of metadata (pdb REMARKs, for example) is specific to the indexing of atoms, residues and chains in that file, which are usually changed by prepmd.

## Running MD simulations with runmd
Steps in the `runmd` workflow:
* Validate user input - runmd will attempt to infer the best parameters and halt if incompatible/impossible settings are used.
* Create an OpenMM system object. If small molecules are present, `runmd` will also load the [OpenFF](https://openforcefield.org/) Sage small molecule force field.
* If there is explicit solvent, set up the simulation box and solvate the system.
* If the run is a metadynamics run, setup bias variables and forces using [openmmtools](https://openmmtools.readthedocs.io/en/stable/).
* Attempt to minimise and run the simulation with OpenMM. If the run/minimisation crashes, the numerical integrator will automatically be switched to the variable langevin integrator and the simulation will be restarted.
* If the run is a metadynamics run, and the metadynamics collective variables aren't minimised, the simulation will restart.

## runmd command-line reference
* Use `runmd --help` for a full list of parameters. 

### Output files
`runmd 6xov_processed.cif -o structure_minimised.cif --traj_out traj.xtc --md_steps 5000` writes out a trajectory file and a structure file (mmcif or pdb) containing the minimised system. If the system has been solvated, this structure file also contains the solvent molecules. The trajectory can be written in DCD or XTC format, which is detected from the filename. The xtc format results in smaller files but with less precision.
### Explicit solvent
`runmd structure.cif -o structure_minimised.cif --traj_out traj.xtc --md_steps 500 --step 10 -solv tip4pew` will run a simulation with the tip4pew solvent. tip3p, tip4pew and spce are supported. By default, simulations run with an implicit solvent equivalent to AMBER's `igb=8` option.
### Temperature and pressure
The default settings result in a rather loose coupling to the heat bath. You can change this with the `-f` or `--friction` argument, which specified the friction coefficient coupling the system to the heat bath. Running a simulation with explicit solvent will result in tighter coupling. You can also add pressure coupling with `--pressure 1.0` (for 1 bar).
### Change force field
`runmd structure.cif -o structure_minimised.cif --traj_out traj.xtc --md_steps 500 --step 50 -ff amber14` runs with amber14. charmm36, amoeba, amber14 and amber19 are available, with charmm36 being the default. The force field is one of the most important MD parameters, and the best force field to use is normally system-dependent.
### Equilibrate only side chains:
`runmd structure.cif -o structure_minimised.cif --fix_backbone -solv tip4pew --notest` will fix the backbone in place and only equilibrate side chains.
### Add ligands
`runmd structure.cif -l LIG.sdf -ff amber14` runs a simulation with a ligand. You can add multiple ligands by using the `-l` argument multiple times. `runmd` supports small molecules using openff's Sage force field, which has limited compatibility with other force fields and solvent models, so ligand simulations only run with the amber14 force field and explicit solvent. By default, ligand simulations also run with a smaller timestep.
### Create a morph trajectory
`runmd pre.cif -m post.cif -o minimised_out.pdb`  will create a trajectory that smoothly transitions between the structures in pre.cif and post.cif. This trajectory is created using openmmtools' metadynamics features. The metadynamics run applies arbitrary biasing forces to perform the transition, so this should only be used for visualisation/illustration, and may not represent the underlying physics and biology.
If you have two files for the same structure which aren't aligned (e.g. they have slightly different starting/ending residues), you can trim the ends to align them: `aligntogether pre.cif post.cif pre_cropped.cif post_cropped.cif`
### Numerical integrators
Set the numerical integrator with the `-i` flag. This can be either `VariableLangevinIntegrator` or `LangevinMiddleIntegrator`. By default, `runmd` will attempt to use the latter, and fall back to the former if the simulation becomes numerically unstable. The parameter `--minimise-err` sets the error tolerance or the variable langevin integrator. Its value is arbitrary - 0.001 is a good starting point, increasing it will make the simulation run faster at the expense of accuracy.
### Other settings
By default, `runmd` will try to select the most optimal nonbonded interaction method, but this can be overridden with `-nb` or `--nonbonded`, which can be one of `PME`, `CutoffPeriodic`, or `CutoffNonPeriodic`. Similarly, it will constrain the length of all bonds involving a hydrogen atom, which can allow for longer timesteps at the cost of some accuracy. This can be disabled by setting `-c None` or `--constraints None`. This setting is also disabled if the backbone is fixed.

## What next?
* Though you can run simple MD simulations, minimisations and validation with `prepmd`, for more in-depth MD we recommend using software such as GROMACS, AMBER, NAMD and OpenMM.
* If you're looking to generate an atomistic structure file that matches your EM map as closely as possible, you can use a flexible fitting tool such as [TEMPy-ReFF](https://gitlab.com/topf-lab/tempy-reff).

## Python API
prepmd is also accessible via a python API:
```python
from prepmd.add_modeller_license import setup_license
setup_license(<your license key>) # only needs to be done once

from prepmd.prep import prep
prep("6xov", "6xov.cif", "working_dir")

from prepmd.run import run
run("6xov.cif", traj_out="traj.xtc")
```

## Licence
AGPLv3

## Contributors
`prepmd` is developed by Rob Welch. Thanks to Harry Swift for helping set up the CI. This project is funded by [DRI-IMB](https://driimb.org/). The repo is managed by CCPBioSim.

## Dependencies
* OpenMM
* PDBFixer
* BioPython
* MODELLER
* pdb2pqr
* mrcfile
* icp
* mdanalysis
* openmmtools
* openff-toolkit
* rdkit

## Developer notes

### Manual install
* Install [Conda](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install) (if you don't already have it)
* Clone this repo and enter the folder: `git clone https://github.com/CCPBioSim/prepmd.git && cd prepmd` 
* Run `conda env create --name prepmd --file environment.yaml && conda activate prepmd && pip install -e .`
* For the MODELLER part of the workflow to work, you need a [modeller license key](https://salilab.org/modeller/registration.html) and add it to modeller's config.py file. If you use conda, the key will be in `envs/prep/lib/modeller-10.7/modlib/modeller/config.py` relative to the path where conda is installed.
* After installing, run `pytest` to run tests.

### Code structure
* Prep structure
    * prep.py - entry point, actual structure preparation happens in prep.prep
    * model.py - calls to the MODELLER API and functions for parsing FASTA files and MODELLER output
    * fix.py - calls to pdbfixer API and other small fixes for structure files
    * get_residues.py - parse mmcif/pdb files and convert them to fasta sequences
    * util.py - utility functions, mostly to do with residue sequences
    * point_cloud.py - compare structure files with EM density maps
    * download.py - download structures from various online services
    * pdb2pqr - call to pdb2pqr API
* Run structure
    * run.py - calls to OpenMM for running simulations
    * metadynamics.py - calls to openmmtools to set up metadynamics simulation
* Other stuff\shared stuff
    * add_modeller_license.py - command-line utility to add a license key to MODELLER
    * ligand.py - parse structures contianing ligands (for system prep) and prepare ligand force field (for running simulations)
    * align_together.py - command-line utility to trim two structure files to be the same length and have corresponding atom indices. Useful
* Tests
    * In test/test_all.py. All integration tests.

## Style
* All code is PEP8 formatted
* Pure functions, minimal coupling and imperative style preferred. Code is mostly WET and functions are allowed to be long to reduce indirection and avoid hiding complexity.
* All user-facing code (e.g. CLI) should have rigorous input validation and descriptive error messages
* prep.py and run.py are both stateful in evil ways; MODELLER writes tons of junk to the working directory, so prep.py changes the working directory to a temporary folder and then copies output files to a path that the user has specified. This necessitates a bit of global state with the working directory, and also juggling different paths and accounting for relative/absolute paths provided by the user.
* run.py is even more stateful - OpenMM objects are pretty stateful, and runmd can re-initialise openmm system objects, keep and restore old co-ordinates (if the simulation becomes numerically unstable) and also recursively call itself to try and fix metadynamics problems. This looks ugly as sin, but it's worth it for those features.