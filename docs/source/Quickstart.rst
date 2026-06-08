Quickstart
==========

`prepmd 6xov 6xov_processed.cif` will download the structure for PDB entry `6xov`, process it and write it to `6xov_processed.cif`.

`runmd 6xov_processed.cif --traj_out traj.xtc --md_steps 5000` will minimise and run a simulation of 6xov_processed.cif writing a trajectory to `traj_out.xtc`, for 5000 steps. By default, `runmd` uses a minimal set of simulation parameters, which aren't likely to be accurate - check the `runmd` section of this documentation for more options.