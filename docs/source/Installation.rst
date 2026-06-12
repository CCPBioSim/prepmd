Installation
============

* Install [Conda](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install) (if you don't already have it). If you have an existing conda install, make sure it can install packages from conda-forge.
* Recommended: create a new virtual environment: `conda env create --name prepmd && conda activate prepmd`
* Install prepmd from the CCPBioSim conda channel: `conda install -c CCPBioSim prepmd`
* Optional: add your [modeller license key](https://salilab.org/modeller/registration.html) by running `prep-license <your license key>`. If you don't have a MODELLER license key, prepmd will instead use pdbfixer's loop builder.