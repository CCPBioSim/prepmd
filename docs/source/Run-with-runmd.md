Running MD simulations with runmd
=================================

## runmd workflow
* Validate user input - runmd will attempt to infer the best parameters and halt if incompatible/impossible settings are used.
* Create an OpenMM system object. If small molecules are present, `runmd` will also load the [OpenFF](https://openforcefield.org/) Sage small molecule force field.
* If there is explicit solvent, set up the simulation box and solvate the system.
* If the run is a metadynamics run, setup bias variables and forces using [openmmtools](https://openmmtools.readthedocs.io/en/stable/).
* Attempt to minimise and run the simulation with OpenMM. If the run/minimisation crashes, the numerical integrator will automatically be switched to the variable langevin integrator and the simulation will be restarted.
* If the run is a metadynamics run, and the metadynamics collective variables aren't minimised, the simulation will restart.

## runmd command-line reference
* Use `runmd --help` for a full list of parameters. 

## Example usage

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
`runmd structure.cif -l LIG.sdf -ff amber14` runs a simulation with a ligand. You can add multiple ligands by using the `-l` argument multiple times. `runmd` supports small molecules using openff's Sage force field, which has limited compatibility with other force fields and solvent models, so ligand simulations only run with the amber14 force field and explicit solvent. By default, ligand simulations also run with a much smaller timestep. Note: runmd will try to infer the ligand name, but you can pass `-ln` for each ligand to specify the name manually. This string will be used to identify the ligand in the output topology file, so if it can't be found, the ligand will just be listed as 'UNK'.

### Create a morph trajectory
`runmd pre.cif -m post.cif -o minimised_out.pdb`  will create a trajectory that smoothly transitions between the structures in pre.cif and post.cif. This trajectory is created using openmmtools' metadynamics features. The metadynamics run applies arbitrary biasing forces to perform the transition, so this should only be used for visualisation/illustration, and may not represent the underlying physics and biology.
If you have two files for the same structure which aren't aligned (e.g. they have slightly different starting/ending residues), you can trim the ends to align them: `aligntogether pre.cif post.cif pre_cropped.cif post_cropped.cif`

### Numerical integrators
Set the numerical integrator with the `-i` flag. This can be either `VariableLangevinIntegrator` or `LangevinMiddleIntegrator`. By default, `runmd` will attempt to use the latter, and fall back to the former if the simulation becomes numerically unstable. The parameter `--minimise-err` sets the error tolerance or the variable langevin integrator. Its value is arbitrary * 0.001 is a good starting point, increasing it will make the simulation run faster at the expense of accuracy.

### Other settings
By default, `runmd` will try to select the most optimal nonbonded interaction method, but this can be overridden with `-nb` or `--nonbonded`, which can be one of `PME`, `CutoffPeriodic`, or `CutoffNonPeriodic`. Similarly, it will constrain the length of all bonds involving a hydrogen atom, which can allow for longer timesteps at the cost of some accuracy. This can be disabled by setting `-c None` or `--constraints None`. This setting is also disabled if the backbone is fixed.
