Developer notes
===============

Dependencies
------------
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

Manual install
--------------
* Install [Conda](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install) (if you don't already have it)
* Clone this repo and enter the folder: `git clone https://github.com/CCPBioSim/prepmd.git && cd prepmd` 
* Run `conda env create --name prepmd --file environment.yaml && conda activate prepmd && pip install -e .`
* For the MODELLER part of the workflow to work, you need a [modeller license key](https://salilab.org/modeller/registration.html) and add it to modeller's config.py file. If you use conda, the key will be in `envs/prep/lib/modeller-10.7/modlib/modeller/config.py` relative to the path where conda is installed.
* After installing, run `pytest` to run tests.

Code structure
--------------
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

Style
-----
* All code is PEP8 formatted
* All user-facing code (e.g. CLI) should have rigorous input validation and descriptive error messages
* Pure functions, minimal coupling and imperative style preferred. Code is mostly WET to reduce coupling where possible, with abstractions only being made where it's obviously necessary to do so. Functions are allowed to be long to reduce indirection and avoid hiding complexity. In particular, longer functions aren't broken up in a way that would require shorter functions to accept objects with a specific state as arguments.
* prep.py and run.py are both stateful in evil ways; MODELLER writes tons of junk to the working directory, so prep.py changes the working directory to a temporary folder and then copies output files to a path that the user has specified. This necessitates a bit of global state, and also juggling different paths and accounting for relative/absolute paths provided by the user.
* run.py is even more stateful - OpenMM objects are pretty stateful, and runmd can re-initialise openmm system objects, keep and restore old co-ordinates (if the simulation becomes numerically unstable) and also recursively call itself to try and fix metadynamics problems. This looks ugly as sin, but it's worth it for those features.