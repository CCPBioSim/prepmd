<p align="center" width="100%">
  <img src="https://raw.githubusercontent.com/CCPBioSim/branding/refs/heads/main/logos/prepmd/prepmd-logo-black-text.svg" width="256px" >
</div>


| Category       | Badges |
|----------------|--------|
| **Build**      | [![prepmd CI](https://github.com/CCPBioSim/prepmd/actions/workflows/python-app.yml/badge.svg)](https://github.com/CCPBioSim/prepmd/actions/workflows/python-app.yml) |
| **Documentation** |[![Docs - Status](https://app.readthedocs.org/projects/prepmd/badge/?version=latest)](https://prepmd.readthedocs.io/en/latest/?badge=latest) |
| **Citation**      |  |
| **Anaconda**       | [![Anaconda.org](https://anaconda.org/CCPBioSim/prepmd/badges/version.svg)](https://anaconda.org/CCPBioSim/prepmd/) [![Last Updated](https://anaconda.org/CCPBioSim/prepmd/badges/latest_release_date.svg)](https://anaconda.org/CCPBioSim/prepmd) [![Platforms](https://anaconda.org/CCPBioSim/prepmd/badges/platforms.svg)](https://anaconda.org/CCPBioSim/prepmd) [![License](https://anaconda.org/CCPBioSim/prepmd/badges/license.svg)](https://anaconda.org/CCPBioSim/prepmd) [![Downloads](https://anaconda.org/CCPBioSim/prepmd/badges/downloads.svg)](https://anaconda.org/CCPBioSim/prepmd)|
| **Quality**    | [![Coverage Status](https://coveralls.io/repos/github/CCPBioSim/prepmd/badge.svg?branch=main)](https://coveralls.io/github/CCPBioSim/prepmd?branch=main) |

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

## Example usage
<img width="768" height="432" alt="demo" src="https://github.com/user-attachments/assets/f7785a4d-ef54-428d-b193-1cb6413bc21e" />

For more detailed info, check the [readthedocs page](https://prepmd.readthedocs.io/en/latest).

## Licence
AGPLv3

## Contributors
`prepmd` is developed by Rob Welch. Thanks to Harry Swift for helping set up the CI, organising a lot of the repo stuff and managing the CCPBioSim conda channel. The logo was created by Jas Kalayan. This project is funded by [DRI-IMB](https://driimb.org/) and the repo is managed by CCPBioSim.
