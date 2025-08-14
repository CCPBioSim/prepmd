# prepmd
A utility to automatically prepare structures from the PDB for molecular dynamics simulation.

## Features
* [X] Automatically download structures, sequences and metadata from the PDB and UNIPROT
* [X] Automatically fill missing loops with modeller
* [X] Automatically add missing atoms and convert non-standard residues (using pdbfixer)
* [X] Automatically propagate metadata through to finalised structure files
* [X] Automatically resolve steric clashes
* [ ] Automatically dock ligands
* [ ] AIIDA integration

## Installation
* Download this repo and enter the root folder
* Run `conda env create --name prepmd --file environment.yaml && conda activate prepmd && pip install .`
* For the modeller part of the workflow to work, you need to get a [modeller license key](https://salilab.org/modeller/registration.html) and add it to modeller's config.py file. If you use conda, the key will be in envs/prepmd/lib/modeller-10.7/modlib/modeller/config.py relative to your conda install directory.
* Run `prempmd` from the prepmd environment. A basic example: `prepmd 6xov 6xov_processed.pdb` will download the structure for PDB entry 6xov, process it and write it to `6xov_processed.pdb`.
