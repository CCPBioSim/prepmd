Preparing structure files for simulation with prepmd
====================================================

Note: When running `prepmd`, we recommend using .mmCif files where possible. The .pdb file format is deprecated and is provided for legacy compatibility only.

## prepmd workflow
Steps in the `prepmd` workflow:
* The structure file(s) are downloaded (if not supplied) into a working directory. `prepmd` automatically infers the file format from the file extension of the input/output files.
* `prepmd` extracts the sequence from the residues in the PDB directly and compares them to a reference sequence. By default this is the sequence described in the SEQRES entries of the structure file. The two sequences are alligned and [MODELLER](https://salilab.org/modeller/) is used to fill in the missing residues.
* Optionally, multiple models can be created, and scored based on MODELLER's internal metrics or their similarity to a reference EM density map.
* Depending on settings, HETATMS may be extracted from the structure file and saved to .sdf files. [rdkit](https://www.rdkit.org/) is used to add hydrogens and correct the geometry of the ligands.
* [PDBFixer](https://github.com/openmm/pdbfixer) is used to add missing hydrogens, remove nonstandard residues and add missing residues, if MODELLER isn't being used.
* Optionally, at this point, a PQR file can be output using [PDB2PQR](https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/apbs/pdb2pqr.html).
* Finally, [OpenMM](https://openmm.org/) is used to perform a test minimisation and simulation. This step ensures that the resulting file is ready for simulation and that there are no steric clashes. If the minimisation or test simulation fails, it will be retried with OpenMM's variable langevin integrator. In testing, this has successfully minimised structure files with high clash scores.
* The final, mimimised structure file will be written out. Note: if ligands are present, the non-minimised structure will be written instead - this is to allow the user to choose which ligand files to include in their final structure, which can be minimised using `runmd`.

## prepmd command-line reference
* Use `prepmd --help` for a full list of parameters. 

## Worked examples

### Using a local structure file
 `prepmd --structure 6xov_input.pdb 6xov 6xov_processed.pdb`. You still need to supply a PDB code, as some of the file formats used by prepmd require one to be present. The code doesn't have to be a 'real' PDB code, e.g. 'AAAA' will work fine. When using this setting, the input and output files must be in the same format * prepmd doesn't perform implicit conversions!

### Generate multiple structure files
`prepmd 6xov 6xov_processed.pdb -n 5` will generate 5 candidate structures and select the best one as determined by MODELLER's internal metrics. Alternatively, `prepmd 6xov 6xov_processed.pdb -n 5 -em 22281 --contour 0.01` will download EMD-22281, the EMDB entry associated with 6XOV, and score the generated models based on their agreement with the EM density map (using the iterative closest point algorithm). The -em setting can also point to a map file. Note: this won't do anything without a MODELLER license key.

### Use refined structures from PDB-REDO
`prepmd 1cbs 1cbs_processed.pdb --redo` will download a refined structure from PDB-REDO, if it is available. Note: not all PDB entries have corresponding PDB-REDO entries.

### Use your own alignments and sequences to fill missing loops
By default, `prepmd` will read missing residues from the pdb/mmcif SEQRES records, attempt to align the missing residues with the currently present residues, and then build missing loops with MODELLER. You can manually provide an aligned FASTA file containing the the complete and incomplete sequences with `--fasta`.  You can also ask prepmd to get the sequence data from UNIPROT instead, with `--download`, though this is not recommended, as the raw sequence data can be substantially different from the PDB and cause the alignment to fail.

### Handling ligands
By default, `prepmd` strips all HETATMS records from structure files. To keep these files, run prepmd with `--hetatms`, which will strip the HETATMS records from the structure file and write each residue to its own SDF file. The co-ordinates inside the SDF files correspond to the co-ordinates of the ligands in the structure file, so the ligands can be added back into the original structure easily. `prepmd` uses [rdkit](https://www.rdkit.org/) to add hydrogens and correct the geometry of small molecules.

### Adding hydrogens
Use the `--ph` flag to set the pH to use when adding hydrogens in pdbfixer, e.g. `prepmd 6xov 6xov.pdb -ph 7.0`. For more control, you can output a PQR file by running `prepmd 6xov 6xov.pdb -ph 7.0 --pqr 6xov.pqr`. As with ligands, this file won't be minimised, as that would remove the charges and radii added by PDB2PQR, though you this file can still be minimised with runmd after the fact.

### Working directory
By default, `prepmd` will leave intermediate files in a randomly-named temporary directory. You can set the name of this directory: `prepmd --wdir 6xov_temp 6xov 6xov.cif`.

### Other notes
Warning: `prepmd`'s output does not contain any of the metadata found in the original pdb. This is an intentional omission - a lot of metadata (pdb REMARKs, for example) is specific to the indexing of atoms, residues and chains in that file, which are usually changed by prepmd.
