
# prepmd
[![prepmd CI](https://github.com/CCPBioSim/prepmd/actions/workflows/python-app.yml/badge.svg)](https://github.com/CCPBioSim/prepmd/actions/workflows/python-app.yml)

A utility to automatically prepare structures from the PDB for molecular dynamics simulation and perform minimisations and simple MD simulations.

## Features
* Automatically download structures, sequences and metadata from the PDB, PDB-REDO, EMDB and UNIPROT
* Automatically fill missing loops with MODELLER or pdbfixer
* Automatically add missing atoms and fix non-standard residues with pdbfixer
* Automatically resolve steric clashes and minimise structures
* Automatically align and trim together structures to be the same length
* Automatically extract and prepare hetatms\ligands for simulation with rdkit
* Easily run simple MD simulations for testing, validation and minimisation with OpenMM
* Create 'morph' trajectories with metadynamics
* Coming soon: integration with other MD\EM workflows!


## Licence
AGPLv3

## Contributors
`prepmd` is developed by Rob Welch. Thanks to Harry Swift for helping set up the CI. This project is funded by [DRI-IMB](https://driimb.org/). The repo is managed by CCPBioSim.
