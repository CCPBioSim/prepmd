#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Use OpenMM to run simulations
"""

from openmm.app import *
from openmm import *
from openmm.unit import *
from operator import methodcaller
import argparse
import sys
import json
import os
import pathlib
import prepmd.metadynamics
import MDAnalysis as mda
from sys import stdout
from prepmd import ligand
from prepmd import analysis

ff_lookup = {
    "amber14,tip3p": ['amber14-all.xml', 'amber14/tip3p.xml'],
    "amber14,spce": ['amber14-all.xml', 'amber14/spce.xml'],
    "amber14,tip4pew": ['amber14-all.xml', 'amber14/tip4pew.xml'],
    "charmm36,tip3p": ['charmm36.xml', 'charmm36/water.xml'],
    "charmm36,tip4pew": ['charmm36.xml', 'charmm36/tip4pew.xml'],
    "charmm36,spce": ['charmm36.xml', 'charmm36/spce.xml'],
    "amoeba,tip3p": ['amoeba2018.xml'],
    "amoeba,tip4pew": ['amoeba2018.xml'],
    "amoeba,spce": ['amoeba2018.xml'],
    "amber19,tip3p": ['amber19-all.xml' 'amber19/tip3p.xml'],
    "amber19,tip4pew": ['amber19-all.xml' 'amber19/tip4pew.xml'],
    "amber19,spce": ['amber19-all.xml' 'amber19/spce.xml'],
    "amber14,None": ['amber14-all.xml'],
    "charmm36,None": ['charmm36.xml'],
    "amoeba,None": ['amoeba2018.xml'],
    "amber19,None": ['amber19-all.xml'],
    "amber19,implicit": ['amber19-all.xml', 'implicit/gbn2.xml'],
    "amber14,implicit": ['amber14-all.xml', 'implicit/gbn2.xml',],
    "charmm36,implicit": ['charmm36.xml', 'implicit/gbn2.xml'],
    "amoeba,implicit": ['amoeba2018.xml', 'amoeba2009_gk.xml'],  
}


def make_forcefield(ff, solvent=None):
    if ff+","+str(solvent) not in ff_lookup.keys():
        err = "Solvent "+str(solvent)+" not in list of available solvents."
        validff = set([item[0] for item in list(
            map(methodcaller("split", ","), ff_lookup.keys()))])
        validsolvent = set([item[1] for item in list(
            map(methodcaller("split", ","), ff_lookup.keys()))])
        err += "Please pick a valid force field and solvent.\n"
        err += "Valid force fields are "+str(validff)+"\n"
        err += "Valid solvents are "+str(validsolvent)+"\n"
        raise ValueError(err)
    ffargs = ff_lookup[ff+","+str(solvent)]
    return ForceField(*ffargs)


def check_platforms():
    i = 0
    platforms = []
    while True:
        try:
            platforms.append(Platform.getPlatform(i).getName())
            i += 1
        except OpenMMException:
            break
    if "CUDA" not in platforms and "HIP" not in platforms and "OpenCL" not in platforms:
        print("WARNING: no GPU platform found. Simulation will be slow!!!")


def test_sim(pdb, no_output = False):
    if no_output:
        min_out = None
    else:
        min_out = pdb
    run(pdb, minimised_structure_out=min_out, md_steps=None,
        integrator_str="LangevinMiddleIntegrator",
        pressure=None, test_sim_steps=250, no_output = no_output)


def run(pdb,
        minimised_structure_out=None,
        traj_out=None,
        max_minimise_iterations=500,
        minimise_error = None, #default = 0.001,
        test_sim_steps=500,
        md_steps=None,
        md_timestep=None, # default = 0.002*picoseconds,
        forcefield_str=None,  # or amber 14, amber19, amoeba
        integrator_str="LangevinMiddleIntegrator",  # or LangevinMiddleIntegrator
        friction_coeff=1/picosecond,
        temperature=300*kelvin,
        minimise=True,
        test_run=True,
        fix_backbone=False,
        constraints_str="Default",  # Default, None, HBonds
        solvent=None,
        strip_solvent=False,
        ionic_strength=0.0*molar,
        pressure=None,# formerly 1*bar,
        step=1000,
        thermo_out_file="thermo.txt",
        non_bonded_method_str="Default",
        checkpoint_output="checkpoint.dat",
        verbose=True,
        write_params="params.json",
        metadynamics_morph = None, # filename of structure to morph to
        meta_rmsd_threshold_nm = 0.17,
        no_output = False,
        ligands = [],
        ligands_names = None
        ):
    """
    Run an MD simulation from a pdb/mmcif structure created with prepmd.
    Args:
        minimised_structure_out - filename to write the final minimised
        structure to (string)
        traj_out - output trajectory filename (xtc or dcd) if MD runs, a string
        max_minimise_iterations - maximum iterations of the energy minimisation
        algorithm, an int
        minimise_error - error tolerance for variable langevin integrator. The
        value is arbitrary, 0.001 is a good starting point, increasing this
        will make the simulation run faster at the expense of accuracy, a float
        test_sim_steps - how many steps to run of the test simulation. This
        isn't production MD, this is just the simulation that checks that your
        structure doesn't have any steric clashes, an int
        md_steps - number of steps for production md, an int
        md_timestep - md simulation timestep. Multiply by an openmm time unit.
        A float.
        forcefield_str - which MD forcefield to use. Valid forcefield: charmm36,
        amber14, amber19 (if using a recent openmm version), amoeba. A str
        integrator_str - which integrator to use. Valid choices are
        LangevinMiddleIntegrator, VariableLangevinIntegrator. Starting with
        the middle integrator is probably best, as the variable langevin
        integrator will be used automatically if the test simulation crashes.
        A string.
        friction_coeff - the friction coefficient which couples the system to
        the heat bath. divide by an openmm time unit. A float
        temperature - simulation temperature, a float. Multiply by
        an openmm temperature unit.
        minimise - whether to run minimisation. You almost certainly want to
        turn this on, unless you are using an already minimised structure.
        A bool.
        test_run - whether to run a short test MD test run. A bool.
        fix_backbone - whether to fix the backbone in place, (for example, if
        you're resolving the positions of side chains). A bool
        constraints_str - whether to constrain Hbonds or not. Possible values:
        None, HBonds. By default, HBonds will be constrained if the backbone
        is not fixed. A string.
        solvent - solvent to use. Possible values: None (no solvent), tip3p, 
        tip4pew, spce, implicit (gbn model, equivalent to AMBER igb=8). A
        string.
        write_solvent: whether to write solvent atoms to the minimised
        structure file. Solvent will always be written to the trajectory. A
        bool.
        ionic_strength - ionic strength of the solvent. A float. multiply by
        openmm molar.
        pressure - pressure coupling via monte carlo barostat. If None, no
        pressure coupling will be used, otherwise specify a pressure multiplied
        by an openmm pressure unit. A float.
        step - how often to write out traj/thermo information. An int.
        thermo_out_file - name of a file to write thermo information to. A
        string.
        non_bonded_method_str - method for calculating long-range interactions. Can
        be one of: PME, CutoffPeriodic, CutoffNonPeriodic. A string. If this is
        not set, it will be picked automatically.
        checkpoint_output: name of checkpoint file to write to. A str.
        write_params: name of a file to write all the simulation params to. a
        str
        ligands - list of ligand filenames. i think? in sdf format.
        ligands_names - a list of strings of names of ligands. used to output a
        ligand topology. Can run without this but ligands will show up as UNK
    """

    # don't look at this
    locals_copy = copy.copy(locals())
    for key, value in locals_copy.items():
        if type(value) == Quantity:
            locals_copy[key] = value.value_in_unit(value.unit)
    locals_json = json.dumps(locals_copy)

    # input validation, errors, warnings etc
    check_platforms()
    if traj_out:
        if ".xtc" not in traj_out.lower() and ".dcd" not in traj_out.lower():
            raise ValueError("Format of output trajectory must be DCD or XTC")

    if md_steps and traj_out and strip_solvent:
        print("Warning: you are removing solvent molecules from the structure"
              " file, but writing a trajectory, which will contain solvent"
              "molecules. Ideally the trajectory and structure file should "
              "have the same number of atoms.")

    if traj_out and not md_steps:
        raise ValueError("Output trajectory specified, but not a number of "
                         "steps!")
    if md_steps and not traj_out:
        raise ValueError("No output trajectory file specified!")

    if (constraints_str == "HBonds" or constraints_str == HBonds) and fix_backbone:
        raise ValueError("Can't fix the backbone with HBonds constraints! "
                         "Please disable constraints or unfix the backbone!")
    
    # set some defaults
    if forcefield_str == None:
        if not ligands:
            forcefield_str = "charmm36"
        else:
            forcefield_str = "amber14"
            
    if solvent == None:
        if not ligands:
            solvent = "implicit"
        else:
            solvent = "tip3p"
            
    # default timestep for ligand simulations is smaller
    if not md_timestep and not ligands:
        md_timestep = 0.002*picoseconds
    if not md_timestep and ligands:
        md_timestep = 0.00025*picoseconds
    if not minimise_error and not ligands:
        minimise_error = 0.001
    if not minimise_error and ligands:
        minimise_error = 0.001*0.25
    #print("Timestep: "+str(md_timestep))
    #print("Min_err: "+str(minimise_error))

    if forcefield_str == "amoeba" and not solvent:
        if non_bonded_method_str == "Default":
            non_bonded_method = NoCutoff
            print("WARNING: AMOEBA does not support non-bonded cutoffs without"
                  " defining a periodic simulation box. For optimal "
                  "performance, use a solvated box or a different force"
                  " field.")
        else:
            raise ValueError("AMOEBA can't be run with the specified non "
                             "bonded method without a simulation box. Try "
                             "running a solvated simulation instead.")

    if not solvent:
        print("WARNING: No solvent. Minimised structures may be wrong!!!")
        
    if fix_backbone and traj_out:
        print("Warning: you are writing a trajectory but have constrained "
              "the backbone. Your system won't move!")
    
    if minimise and not minimised_structure_out and not no_output:
        minimised_structure_out = "min_"+os.path.basename(pdb)
        print("No minimised structure output file specified. Writing "
              "minimised structure to "+"min_"+pdb)

    if ligands:
        if forcefield_str != "amber14":
            raise ValueError("Error: Ligand simulations are only supported "
                             "with amber14")
        if solvent != "tip3p" and solvent != "tip4pew":
            raise ValueError("Ligand simulations are only supported with"
                             " tip3p and tip4pew solvent")
        if strip_solvent:
            raise ValueError("Can't strip solvent on a ligand system!")
        if constraints_str == "Default" or constraints_str == "HBonds":
            ff = "ff14sb_off_impropers_0.0.4.offxml"
            ligand_ff = "openff-2.2.0.offxml"
        if constraints_str == "None":
            ff = "ff14sb_off_impropers_0.0.4.offxml"
            ligand_ff = "openff_unconstrained-2.2.0.offxml"
        if md_timestep >= 0.001*picoseconds or minimise_error >= 0.001:
            print("WARNING: Timestep/error tolerance are high for a ligand "
                  "simulation. If the simulation crashes (e.g. due to "
                  "'Particle co-ordinate is NaN'), try reducing them.")
        if ".cif" in pdb or ".mmcif" in pdb:
                writer = PDBxFile
        else:
                writer = PDBFile
        system, topology, positions = ligand.create_ligand_system(pdb,
                                                                  ligands,
                                                                  water_ff = solvent,
                                                                  desired_charge_elementary = ionic_strength._value,
                                                                  ff = ff,
                                                                  ligand_ff = ligand_ff,
                                                                  n_water = 5000,
                                                                  ligands_names=ligands_names)
        print("Note: topology written to "+str(minimised_structure_out)+" is "
              "not minimised. This is because OpenMM cannot write a correct "
              "topology for OpenFF ligands and solvents.")
        # TODO: maybe write a separate file?
    
    else:
        # read in mmcif or pdb
        if ".cif" in pdb or ".mmcif" in pdb:
            structure = PDBxFile(pdb)
            writer = PDBxFile
        else:
            structure = PDBFile(pdb)
            writer = PDBFile
    
        modeller = Modeller(structure.topology, structure.positions)
    
        forcefield = make_forcefield(forcefield_str, solvent)
    
        if solvent and solvent != "implicit":
            print("Solvating system...")
            modeller.addSolvent(forcefield, model=solvent,
                                ionicStrength=ionic_strength,
                                padding=1.0*nanometers)
    
        # can't set hbond constraints if also restraining whole backbone
        if constraints_str not in [None, "Default", HBonds, "HBonds"]:
            raise ValueError(
                "Constraints must be one of: Default, None, HBonds")
        if constraints_str == "Default":
            if fix_backbone:
                constraints = None
                print("As the backbone is fixed, no constraints will be "
                      "applied.")
            else:
                constraints = HBonds
                print("No constraints specified. Constraining HBonds.")
        if constraints_str == "HBonds":
            constraints = HBonds
    
        # setup non bonded method - if none is chosen, select based on box
        if non_bonded_method_str == "Default":
            if solvent == None or solvent == "implicit":
                non_bonded_method = CutoffNonPeriodic
            else:
                non_bonded_method = PME
        elif non_bonded_method_str == "PME":
            non_bonded_method = PME
        elif non_bonded_method_str == "CutoffPeriodic":
            non_bonded_method = CutoffPeriodic
        elif non_bonded_method_str == "CutoffNonPeriodic":
            non_bonded_method = CutoffNonPeriodic
        elif non_bonded_method == NoCutoff:
            pass
        else:
            raise ValueError("Supply a valid nonbondedMethod! Could be "
                             "Default, PME, CutoffPeriodic or "
                             "CutoffNonPeriodic")
        if non_bonded_method == PME and not solvent:
            raise ValueError("Cannot run a PME simulation on a system with no "
                             "simulation box. Either change the "
                             "nonbondedMethod to 'NoCutoff' or add a solvent,"
                             " which will initialise a box automatically.")
        
        #modeller.addHydrogens()
        system = forcefield.createSystem(modeller.topology,
                                         nonbondedMethod=non_bonded_method,
                                         nonbondedCutoff=1*nanometer,
                                         constraints=constraints)

    # constrain backbone by setting CA mass to zero
    if fix_backbone:
        for atom in modeller.topology.atoms():
            if atom.name == 'CA':
                system.setParticleMass(atom.index, 0*amu)

    if pressure:
        print("Adding pressure coupling: "+str(pressure))
        system.addForce(MonteCarloBarostat(pressure, temperature))
        
    if metadynamics_morph:
        final = PDBFile(metadynamics_morph)
        ref_positions = final.getPositions(asNumpy=False)
        ref_positions_np = final.getPositions(asNumpy=True)
        r = RMSDForce(ref_positions)
        rforce = CustomCVForce('1*r')
        rforce.addCollectiveVariable('r', r)
        rmsd = BiasVariable(force=rforce, minValue=0 * angstrom,
                            maxValue=500 * angstrom, biasWidth=1 * angstrom,
                            periodic=False)
        meta = Metadynamics(system=system, variables=[rmsd],
                            temperature=300, biasFactor=2500.0,
                            height=10 * kilocalories_per_mole,
                            frequency=5, saveFrequency=50, biasDir='./')

    # set up integrator
    if integrator_str == "VariableLangevinIntegrator":
        integrator = VariableLangevinIntegrator(temperature,
                                                friction_coeff, minimise_error)
    elif integrator_str == "LangevinMiddleIntegrator":
        integrator = LangevinMiddleIntegrator(temperature,
                                              friction_coeff, md_timestep)
    else:
        raise ValueError("Integrator must be one of: "
                         "VariableLangevinIntegrator,"
                         " LangevinMiddleIntegrator")
        
    if ligands:
        simulation = Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)
    else:
        topology = modeller.topology
        positions = modeller.positions
        simulation = Simulation(topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

    # run simulation
    # try and resolve steric clashes with variable langevin integrator
    # todo: this section could be tidier i think
    strict_min = quantity.Quantity(value=2.5,unit=kilojoule/(nanometer*mole))
    if minimise:
        try:
            print("Minimising...")
            simulation.minimizeEnergy(maxIterations=max_minimise_iterations,
                                      tolerance = strict_min)
            curr_state = simulation.context.getState(
                getPositions=True).getPositions(asNumpy=True)
        except OpenMMException as e:
            if integrator_str != "VariableLangevinIntegrator":
                print("Minimisation blew up. "
                      "Trying variable langevin integrator...")
                integrator = VariableLangevinIntegrator(temperature,
                                                        friction_coeff,
                                                        minimise_error)
                simulation = Simulation(topology, system, integrator)
                simulation.context.setPositions(positions)
                integrator_str = "VariableLangevinIntegrator"
                simulation.minimizeEnergy(
                    maxIterations=max_minimise_iterations,
                    tolerance = strict_min)
                curr_state = simulation.context.getState(
                    getPositions=True).getPositions(asNumpy=True)
            else:
                print("Could not minimise even with variable langevin "
                      "integrator.")
                raise e
    
    if test_run:
        try:
            print("Running test simulation...")
            simulation.step(test_sim_steps)
            print("Test simulation ran successfully.")
        except OpenMMException as e:
            if integrator_str != "VariableLangevinIntegrator":
                print("Simulation blew up, running with variable langevin "
                      "integrator...")
                integrator = VariableLangevinIntegrator(temperature,
                                                        friction_coeff,
                                                        minimise_error)
                simulation = Simulation(topology, system, integrator)
                simulation.context.setPositions(positions)
                simulation.minimizeEnergy(
                    maxIterations = max_minimise_iterations,
                    tolerance = strict_min)
                print("Fixed.")
                curr_state = simulation.context.getState(
                    getPositions=True).getPositions(asNumpy=True)
                if md_steps:
                    print("Note: production MD will run with the variable "
                          "langevin integrator.")
            else:
                print("Simulation blew up even while using variable langevin"
                      " integrator. This doesn't normally happen.")
                raise e

    with open(write_params, "w") as file:
        file.write(locals_json)

    if (minimise or test_run) and not no_output:
        modeller_out = Modeller(simulation.topology, curr_state)
        if strip_solvent:
            modeller_out.deleteWater()
        writer.writeFile(modeller_out.topology, modeller_out.positions,
                         file=open(minimised_structure_out, "w"), keepIds=True)
        print("Wrote minimised structure to "+str(minimised_structure_out)+".")
    else:
        if minimised_structure_out and not (minimise or test_run):
            raise ValueError("Minimised structure output was requested, but"
                             " minimisation and test run were both skipped.")
        
    if md_steps and traj_out and not metadynamics_morph:
        if ".dcd" in traj_out.lower():
            simulation.reporters.append(DCDReporter(traj_out, step))
        elif ".xtc" in traj_out.lower():
            simulation.reporters.append(XTCReporter(traj_out, step))
        if thermo_out_file:
            simulation.reporters.append(StateDataReporter(thermo_out_file,
                                                          step, step=True,
                                                          potentialEnergy=True,
                                                          kineticEnergy=True,
                                                          elapsedTime=True,
                                                          temperature=True))
        if verbose:
            simulation.reporters.append(StateDataReporter(sys.stdout,
                                                          step, step=True,
                                                          speed=True,
                                                          remainingTime=True,
                                                          totalSteps=md_steps,
                                                          potentialEnergy=True,
                                                          temperature=True))

    if md_steps and traj_out and not metadynamics_morph:     
        print("Running production MD...")
        simulation.step(md_steps)
        simulation.saveCheckpoint(checkpoint_output)
        print("Done!")
        print("Wrote trajectory to "+traj_out)
        if thermo_out_file:
            print("Wrote thermo info to "+thermo_out_file)
        print("Wrote checkpoint to "+checkpoint_output)
        
    if ligand and traj_out and minimised_structure_out:
        print("Running ligand analysis...")
        analysis.get_ligand_centroid_traj(minimised_structure_out, traj_out)
        print("Done.")
        
    # metadynamics run
    
    if metadynamics_morph and not traj_out:
        print("Metadynamics run requested but no output trajectory specified. "
              "Writing to metadynamics.xtc.")
        traj_out = "metadynamics.xtc"
    if metadynamics_morph and not md_steps:
        print("Metadynamics run requested but no steps specified. Will run a "
              "maximum of 5000000 steps.")
        md_steps = 5000000
    
    write_morph = False
    redo_run = False
    if metadynamics_morph:
        try:
            # Set up the integrator and simulation reporters
            simulation.reporters.append(XTCReporter(traj_out, 500))
            simulation.reporters.append(StateDataReporter(stdout, 500, step=True,
                                                          potentialEnergy=True, temperature=True, progress=True, remainingTime=True, totalSteps=5000000, separator='\t'))
            simulation.reporters.append(prepmd.metadynamics.RMSDReporter(
                "rmsd.txt", 500, ref_positions_np, threshold_nm=meta_rmsd_threshold_nm))
            meta.step(simulation, md_steps)
            print("WARNING: Couldn't reach target structure.")
        except prepmd.metadynamics.StopSimulation:
            print("Metadynamics run finished. Writing morph data...")
            write_morph = True
            pass
        except prepmd.metadynamics.RMSDIncrease:
            print("RMSD is increasing, restarting...")
            redo_run = True
            pass

        if redo_run:
            run(pdb, minimised_structure_out, traj_out,
                max_minimise_iterations, minimise_error,
                test_sim_steps, md_steps, md_timestep, forcefield_str,
                integrator_str, friction_coeff, temperature, minimise,
                test_run, fix_backbone, constraints_str, solvent,
                strip_solvent, ionic_strength, pressure, step,
                thermo_out_file, non_bonded_method_str, checkpoint_output,
                verbose, write_params, metadynamics_morph)
    
        if write_morph:
            # write out pdbs and cifs as well
            final_traj = mda.Universe(pdb, traj_out)
            frames = prepmd.metadynamics.get_representative_rmsd_frames(
                simulation.reporters[2].rmsd_log, 10)
            filtered = final_traj.trajectory[frames]
            with mda.Writer("morph.xtc", final_traj.atoms.n_atoms) as w:
                for frame in filtered:
                    w.write(final_traj)
            for frame in frames:
                with mda.Writer("snapshot-"+str(frame)+".pdb") as w:
                    final_traj.trajectory[frame]
                    w.write(final_traj)
            print("Done simulating.")
        


def entry_point():
    "CLI entry point function. Uses sys.argv and argparse args object."

    parser = argparse.ArgumentParser(prog="runmd",
                                     description="Run MD simulations")
    parser.add_argument("pdb", help="Input filename (pdb or mmcif)")
    parser.add_argument(
        "-o", "--min_out", help="filename to write the final minimised structure to (pdb or mmCif)", default=None)
    parser.add_argument(
        "-tr", "--traj_out", help="Production MD trajectory output file. Can be in DCD or XTC format.", default=None)
    parser.add_argument("-it", "--minimise_iterations",
                        help="Maximum minimise iterations", default=100, type=int)
    parser.add_argument("-e", "--minimise_err", help="error tolerance for variable langevin integrator. The value is arbitrary, 0.001 is a good starting point, increasing this will make the simulation run faster at the expense of accuracy", default=None, type=float)
    parser.add_argument("-test", "--test_steps", help="how many steps to run of the test simulation. This isn't production MD, this is just the simulation that checks that your structure doesn't have any steric clashes", default=500, type=int)
    parser.add_argument("-mdsteps", "--md_steps",
                        help="how many steps to run production MD for", default=None, type=int)
    parser.add_argument(
        "-ts", "--timestep", help="MD timestep in picoseconds", default=None, type=float)
    parser.add_argument("-ff", "--forcefield",
                        help="Force field to use (can be charmm36, amber14, amoeba, or amber19 on newer openmm versions", default=None)
    parser.add_argument("-i", "--integrator", help="Numerical integration scheme. Can be one of: VariableLangevinIntegrator, LangevinMiddleIntegrator",
                        default="LangevinMiddleIntegrator")
    parser.add_argument(
        "-f", "--friction", help="Friction coefficient which couples the system to the heat bath. Units 1/picosecond.", default=1, type=float)
    parser.add_argument("-t", "--temperature",
                        help="Simulation temperature in Kelvin.", default=300, type=float)
    parser.add_argument("-nomin", "--no_minimise",
                        help="Skip minimisation", action="store_true")
    parser.add_argument("-notest", "--no_testrun",
                        help="Skip test run", action="store_true")
    parser.add_argument("-fb", "--fix_backbone",
                        help="Fix the backbone", action="store_true")
    parser.add_argument("-c", "--constraints", help="Constraints to apply. Possible values: None, HBonds. Constraining HBonds can make simulations run faster, at the cost of some accuracy. You can't constrain HBonds if you're also fixing the backbone, for obvious reasons.", default="Default")
    parser.add_argument(
        "-solv", "--solvent", help="Solvent to use. Possible values: tip3p, tip4pew, spce, implicit (GBn model, equivalent to AMBER igb=8). If no solvent is specified, implicit will be used.", default=None)
    parser.add_argument("-ns", "--strip_solvent",
                        help="Don't write solvent molecules to the minimised output file", action="store_true")
    parser.add_argument("-ion", "--ionic_strength",
                        help="Ionic strength of the solvent in molar.", default=0.0, type=float)
    parser.add_argument(
        "-p", "--pressure", help="Pressure in bar, coupling via monte carlo barostat If pressure is not specified, no pressure coupling will be used.", default=None)
    parser.add_argument(
        "-step", "--step", help="In production MD, write to output files (traj, thermo) every x steps.", default=1000, type=int)
    parser.add_argument(
        "-th", "--thermo", help="Thermo information output file.", default="thermo.txt")
    parser.add_argument(
        "-nb", "--nonbonded", help="Non-bonded interaction method. Can be: PME, CutoffPeriodic, CutoffNonPeriodic.", default="Default")
    parser.add_argument("-ckpt", "--checkpoint",
                        help="Checkpoint output filename", default="checkpoint.dat")
    parser.add_argument(
        "-q", "--quiet", help="Don't print info to the stdout", action="store_true")
    parser.add_argument("-parm", "--write_params",
                        help="Write simulation params to a (json) file", default="params.json")
    parser.add_argument("-m", "--metamorph",
                        help="Instead of running full MD, create a morph trajectory", default=None)
    parser.add_argument("-l", "--ligand", action="append",
                        help="Add a ligand from an sdf file. Call multiple times to add multiple ligands.", default=None)
    parser.add_argument("-ln", "--ligname", action="append",
                        help="Name the ligand. This is how ligands will be named in the topology file. Call once for each ligand.", default=None)

    args = parser.parse_args()
    minimise = not args.no_minimise
    testrun = not args.no_testrun
    verbose = not args.quiet
    run(args.pdb,
        minimised_structure_out=args.min_out,
        traj_out=args.traj_out,
        max_minimise_iterations=args.minimise_iterations,
        minimise_error=args.minimise_err,
        test_sim_steps=args.test_steps,
        md_steps=args.md_steps,
        md_timestep=args.timestep*picoseconds,
        forcefield_str=args.forcefield,
        integrator_str=args.integrator,
        friction_coeff=args.friction/picosecond,
        temperature=args.temperature*kelvin,
        minimise=minimise,
        test_run=testrun,
        fix_backbone=args.fix_backbone,
        constraints_str=args.constraints,
        solvent=args.solvent,
        strip_solvent=args.strip_solvent,
        ionic_strength=args.ionic_strength*molar,
        pressure=args.pressure*bar,
        step=args.step,
        thermo_out_file=args.thermo,
        non_bonded_method_str=args.nonbonded,
        checkpoint_output=args.checkpoint,
        verbose=verbose,
        write_params=args.write_params,
        metadynamics_morph=args.metamorph,
        ligands=args.ligand,
        ligands_names=args.ligname
        )


if __name__ == "__main__":
    entry_point()
