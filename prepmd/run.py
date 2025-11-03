#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Use openmm to run simulations
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
#from sys import stdout
#import time
#import json
#import sys

#forcefields = {
#    "charmm36": ForceField('charmm36.xml', 'charmm36/water.xml'),
#    "amber14": ForceField('amber14-all.xml', 'amber14/tip4pew.xml', 'amber14/spce.xml'),
#    "amoeba": ForceField('amoeba2018.xml'),
#    #"amber19": ForceField('amber19-all.xml', 'amber19/tip3pfb.xml')
#}

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
    "amber19,None": ['amber19-all.xml']
}

def make_forcefield(ff, solvent=None):
    if ff+","+str(solvent) not in ff_lookup.keys():
        validff = set([item[0] for item in list(
            map(methodcaller("split", ","), ff_lookup.keys()))])
        validsolvent = set([item[1] for item in list(
            map(methodcaller("split", ","), ff_lookup.keys()))])
        err = "Please pick a valid force field and solvent.\n"
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


# sarah says to use the ABC protocols - i still don't know what those are

# def run_nounits(pdb,
#         minimised_structure_out=None,
#         traj_out=None,
#         max_minimise_iterations=50,
#         minimise_error=0.001,
#         test_sim_steps=500,
#         md_steps=10000000,
#         md_timestep_picoseconds=0.002,
#         forcefield="charmm36",
#         integrator="LangevinMiddleIntegrator",
#         friction_coeff_per_picosecond=1,
#         temperature_kelvin=300,
#         minimise=True,
#         test_run=True,
#         fix_backbone=False,
#         constraints="Default",
#         solvent="tip3p",
#         write_solvent=False,
#         ionic_strength_molar=0.1,
#         pressure_bar=1,
#         step=1000,
#         thermo_out_file="thermo.txt",
#         non_bonded_method = "Default",
#         checkpoint_output = "checkpoint.dat"
#     ):
    
#     run(pdb,
#         minimised_structure_out,
#         traj_out,
#         max_minimise_iterations,
#         minimise_error,
#         test_sim_steps,
#         md_steps,
#         md_timestep=md_timestep_picoseconds*picoseconds,
#         forcefield=forcefield,
#         integrator=integrator,
#         friction_coeff=friction_coeff_per_picosecond/picosecond,
#         temperature=temperature_kelvin*kelvin,
#         minimise=minimise,
#         test_run=test_run,
#         fix_backbone=fix_backbone,
#         constraints=constraints,
#         solvent=solvent,
#         write_solvent=write_solvent,
#         ionic_strength=ionic_strength_molar*molar,
#         pressure=pressure_bar*bar,
#         step=step,
#         thermo_out_file=thermo_out_file,
#         non_bonded_method = non_bonded_method,
#         checkpoint_output = checkpoint_output
#     )


def test_sim(pdb):
    run(pdb, minimised_structure_out=pdb, md_steps = None,
        integrator = "LangevinMiddleIntegrator", solvent = None,
        pressure = None, test_sim_steps=50)
    
 # todo: why does it not minimise/testrun by default?
 # and why is the backbone constrained by default?

def run(pdb,
        minimised_structure_out=None,
        traj_out=None,
        max_minimise_iterations=100,
        minimise_error=0.001,
        test_sim_steps=500,
        md_steps=None,
        md_timestep=0.002*picoseconds,
        forcefield="charmm36", # or amber 14, amber19, amoeba
        integrator="LangevinMiddleIntegrator", # or LangevinMiddleIntegrator
        friction_coeff=1/picosecond,
        temperature=300*kelvin,
        minimise=True,
        test_run=True,
        fix_backbone=False,
        constraints="Default", #Default, None, HBonds
        solvent="tip4pew",
        strip_solvent=False,
        ionic_strength=0.1*molar,
        pressure=1*bar,
        step=1000,
        thermo_out_file="thermo.txt",
        non_bonded_method = "Default",
        checkpoint_output = "checkpoint.dat",
        verbose = True,
        write_params = "params.json"
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
        forcefield - which MD forcefield to use. Valid forcefield: charmm36,
        amber14, amber19 (if using a recent openmm version), amoeba. A str
        integrator - which integrator to use. Valid choices are
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
        constraints - whether to constrain Hbonds or not. Possible values:
        None, HBonds. By default, HBonds will be constrained if the backbone
        is not fixed. A string.
        solvent - solvent to use. Possible values: None (no solvent), tip3p, 
        tip4pew, spce. A string.
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
        non_bonded_method - method for calculating long-range interactions. Can
        be one of: PME, CutoffPeriodic, CutoffNonPeriodic. A string. If this is
        not set, it will be picked automatically.
        checkpoint_output: name of checkpoint file to write to. A str.
        write_params: name of a file to write all the simulation params to. a
        str
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
        
    if (constraints == "HBonds" or constraints == HBonds) and fix_backbone:
        raise ValueError("Can't fix the backbone with HBonds constraints! "
                         "Please disable constraints or unfix the backbone!")
        
    if forcefield == "amoeba" and not solvent:
        if non_bonded_method == "Default":
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
        print("WARNING: No solvent")

    # read in mmcif or pdb
    if ".cif" in pdb or ".mmcif" in pdb:
        structure = PDBxFile(pdb)
        writer = PDBxFile
    else:
        structure = PDBFile(pdb)
        writer = PDBFile
        
    modeller = Modeller(structure.topology, structure.positions)
    
    forcefield = make_forcefield(forcefield, solvent)
    
    if solvent:
        print("Solvating system...")
        modeller.addSolvent(forcefield, model=solvent,
                            ionicStrength=ionic_strength,
                            padding=1.0*nanometers)
    
    # can't set hbond constraints if also restraining whole backbone
    if constraints not in [None, "Default", HBonds, "HBonds"]:
        raise ValueError(
                         "Constraints must be one of: Default, None, HBonds")
    if constraints == "Default":
        if fix_backbone:
            constraints = None
            print("As the backbone is fixed, no constraints will be applied.")
        else:
            constraints = HBonds
            print("No constraints specified. Constraining HBonds.")
    if constraints == "HBonds":
        constraints = HBonds
    
    # setup non bonded method - if none is chosen, select based on box
    if non_bonded_method == "Default":
        if solvent is None:
            non_bonded_method = CutoffNonPeriodic
        else:
            non_bonded_method = PME
    elif non_bonded_method == "PME":
        non_bonded_method = PME
    elif non_bonded_method == "CutoffPeriodic":
        non_bonded_method = CutoffPeriodic
    elif non_bonded_method == "CutoffNonPeriodic":
        non_bonded_method = CutoffNonPeriodic
    elif non_bonded_method == NoCutoff:
        pass
    else:
        raise ValueError("Supply a valid nonbondedMethod! Could be Default, "
                         "PME, CutoffPeriodic or CutoffNonPeriodic")
    if non_bonded_method == PME and not solvent:
        raise ValueError("Cannot run a PME simulation on a system with no "
                         "simulation box. Either change the nonbondedMethod"
                         "to 'NoCutoff' or add a solvent, which will"
                         "initialise a box automatically.")
        
    system = forcefield.createSystem(modeller.topology,
                                     nonbondedMethod=non_bonded_method,
                                     nonbondedCutoff=1*nanometer,
                                     constraints=constraints)
    
    # constrain backbone by setting CA mass to zero
    if fix_backbone:
        for atom in modeller.topology.atoms():
            if atom.name == 'CA':
                system.setParticleMass(atom.index, 0*amu)
    
    if fix_backbone and traj_out:
        print("Warning: you are writing a trajectory but have constrained "
              "the backbone. Your system won't move!")
    
    if pressure:
        # TODO: this seems to be broken on my current OMM version
        system.addForce(MonteCarloBarostat(pressure, temperature))
    
    # set up integrator
    if integrator == "VariableLangevinIntegrator":
        integrator = VariableLangevinIntegrator(temperature,
                                                friction_coeff, minimise_error)
    elif integrator == "LangevinMiddleIntegrator":
        integrator = LangevinMiddleIntegrator(temperature,
                                              friction_coeff, md_timestep)
    else:
        raise ValueError("Integrator must be one of: "
                         "VariableLangevinIntegrator,"
                         " LangevinMiddleIntegrator")
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    
    # run simulation
    # try and resolve steric clashes with variable langevin integrator
    if minimise:
        print("Minimising...")
        simulation.minimizeEnergy(maxIterations=max_minimise_iterations)
        curr_state = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
    if test_run:
        try:
            print("Running test simulation...")
            simulation.step(test_sim_steps)
            print("Test simulation ran successfully.")
        except OpenMMException as e:
            if integrator != "VariableLangevinIntegrator":
                print("Simulation blew up, running with variable langevin "
                      "integrator...")
                integrator = VariableLangevinIntegrator(temperature,
                                                        friction_coeff,
                                                        minimise_error)
                simulation = Simulation(modeller.topology, system, integrator)
                simulation.context.setPositions(modeller.positions)
                simulation.minimizeEnergy(
                    maxIterations=max_minimise_iterations)
                print("Fixed")
                curr_state = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
                if md_steps:
                    print("Note: production MD will run with the variable"
                          "langevin integrator.")
            else:
                print("Simulation blew up even while using variable langevin"
                      " integrator. This doesn't normally happen.")
                raise e
                
    with open(write_params, "w") as file:
        file.write(locals_json)
    
    if minimise or test_run:
        modeller_out = Modeller(simulation.topology, curr_state)
        if strip_solvent:
            modeller_out.deleteWater()
        writer.writeFile(modeller_out.topology, modeller_out.positions,
            file=open(minimised_structure_out, "w"), keepIds=True)
        print("Wrote minimised structure to "+str(minimised_structure_out)+".")
    else:
        print("Skipped minimisation and test run.")
        if minimised_structure_out:
            raise ValueError("Minimised structure output was requested, but"
                             " minimisation and test run were both skipped.")
    
    if md_steps and traj_out:
        if ".dcd" in traj_out.lower():
            simulation.reporters.append(DCDReporter(traj_out, step))
        elif ".xtc" in traj_out.lower():
            simulation.reporters.append(XTCReporter(traj_out, step))
#        else:
#            raise ValueError("Format of output trajectory must be DCD or XTC")
        if thermo_out_file:
            simulation.reporters.append(StateDataReporter(thermo_out_file,
                                                          step, step=True,
                                                          potentialEnergy=True,
                                                          temperature=True))
        if verbose:
            simulation.reporters.append(StateDataReporter(sys.stdout,
                                        step, step=True,
                                    potentialEnergy=True, temperature=True))
        print("Running production MD...")
        simulation.step(md_steps)
        simulation.saveCheckpoint(checkpoint_output)
        print("Done!")
        print("Wrote trajectory to "+traj_out)
        if thermo_out_file:
            print("Wrote thermo info to "+thermo_out_file)
        print("Wrote checkpoint to "+checkpoint_output)


# no longer used, will be removed soon
# use pytest instead

# def tests():
    
#     file_path = os.path.realpath(__file__)
#     sep = os.path.sep
#     test_file = os.path.dirname(file_path)+sep+"test_data"+sep+"101M_proc.cif"
#     out_dir = "/tmp/prepmdtest/"
#     results = {}
    
#     if not os.path.isdir(out_dir):
#         pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
        
#     print('path =', test_file)
    
    
#     try:
#         run(test_file,
#             traj_out=out_dir+"101M_proc.xtc", md_steps=100, step=25,
#             solvent=None, pressure=None, minimise=False, test_run=False,
#             write_params=out_dir+sep+"params.json",
#             thermo_out_file=out_dir+sep+"thermo.txt",
#             checkpoint_output=out_dir+sep+"checkpoint.dat")
#         results["Basic"] = "Ran successfully"
#     except Exception as e:
#         results["Basic"] = str(e)

#     try:
#         run(test_file, minimised_structure_out=out_dir+"101m_prof_min.cif",
#             traj_out=out_dir+"101M_proc.xtc", md_steps=100, step=25,
#             solvent=None, pressure=None,
#             write_params=out_dir+sep+"params.json",
#             thermo_out_file=out_dir+sep+"thermo.txt",
#             checkpoint_output=out_dir+sep+"checkpoint.dat")
#         results["Minimise/test run"] = "Ran successfully"
#     except Exception as e:
#         results["Minimise/test run"] = str(e)

#     try:
#         run(test_file,
#             traj_out=out_dir+"101M_proc.xtc", md_steps=10, step=2,
#             integrator="VariableLangevinIntegrator",
#             solvent="tip4pew", pressure=1.0*bar, minimise=False,
#             test_run=False,
#             write_params=out_dir+sep+"params.json",
#             thermo_out_file=out_dir+sep+"thermo.txt",
#             checkpoint_output=out_dir+sep+"checkpoint.dat")
#         results["Variable langevin, tip4pew, pressure"] = "Ran successfully"
#     except Exception as e:
#         results["Variable langevin, tip4pew, pressure"] = str(e)

#     try:
#         run(test_file,
#             traj_out=out_dir+"101M_proc.xtc", md_steps=5, step=1,
#             solvent="tip3p", forcefield="amber14", minimise=False,
#             test_run=False,
#             write_params=out_dir+sep+"params.json",
#             thermo_out_file=out_dir+sep+"thermo.txt",
#             checkpoint_output=out_dir+sep+"checkpoint.dat")
#         results["Amber14"] = "Ran successfully"
#     except Exception as e:
#         results["Amber14"] = str(e)
    
#     try:
#         run(test_file,
#             traj_out=out_dir+"101M_proc.xtc", md_steps=5, step=1,
#             solvent="tip4pew", pressure=1.0*bar, forcefield="amber14",
#             fix_backbone=True, minimise=False, test_run=False,
#             write_params=out_dir+sep+"params.json",
#             thermo_out_file=out_dir+sep+"thermo.txt",
#             checkpoint_output=out_dir+sep+"checkpoint.dat")
#         results["Fix backbone"] = "Ran successfully"
#     except Exception as e:
#         results["Fix backbone"] = str(e)
        
#     print("")
#     print("Test results:")
        
#     for name, result in results.items():
#         print(name+": "+result)
    
#    runmd 101M_proc.cif -o 101M_proc_min.cif --traj_out 101M_proc.xtc --md_steps 500 --step 50

#    runmd 101M_proc.cif -o 101M_proc_min.cif --traj_out 101M_proc.xtc --md_steps 5000 --step 50 --nomin --notest # error - minimised structure requested
    
    # run from an already-minimised structure
#    runmd 101M_proc.cif --traj_out 101M_proc.xtc --md_steps 5000 --step 50 -nomin -notest
    
    # solvated
#    runmd 101M_proc.cif -o 101M_proc_min.cif --traj_out 101M_proc.xtc --md_steps 500 --step 50 -solv tip4pew
    
    # solvated with pressure # TODO: seems to be broken on my openmm version (which is a nightly build i think)
#    runmd 101M_proc.cif -o 101M_proc_min.cif --traj_out 101M_proc.xtc --md_steps 500 --step 50 -solv tip4pew -p 1.0
    
    # run with amber14 instead
#    runmd 101M_proc.cif -o 101M_proc_min.cif --traj_out 101M_proc.xtc --md_steps 500 --step 50 -ff amber14
    
    # fix the backbone and minimise to resolve the side chain positions
#    runmd 101M_proc.cif -o 101M_proc_min.cif --fix_backbone -solv tip4pew --notest
    
#    sys.exit()


def entry_point():
    "CLI entry point function. Uses sys.argv and argparse args object."
    #if len(sys.argv) == 2:
    #    if sys.argv[1] == "test":
    #        tests()
            
    parser = argparse.ArgumentParser(prog="runmd",
                        description="Run MD simulations")
                        #epilog="Also, run 'runmd test' for the test suite")
    parser.add_argument("pdb", help="Input filename (pdb or mmcif)")
    parser.add_argument("-o", "--min_out", help="filename to write the final minimised structure to (pdb or mmCif)", default=None)
    parser.add_argument("-tr", "--traj_out", help="Production MD trajectory output file. Can be in DCD or XTC format.", default=None)
    parser.add_argument("-it", "--minimise_iterations", help="Maximum minimise iterations", default=100, type=int)
    parser.add_argument("-e", "--minimise_err", help="error tolerance for variable langevin integrator. The value is arbitrary, 0.001 is a good starting point, increasing this will make the simulation run faster at the expense of accuracy", default=0.001, type=float)
    parser.add_argument("-test", "--test_steps", help="how many steps to run of the test simulation. This isn't production MD, this is just the simulation that checks that your structure doesn't have any steric clashes", default=500, type=int)
    parser.add_argument("-mdsteps", "--md_steps", help="how many steps to run production MD for", default=None, type=int)
    parser.add_argument("-ts", "--timestep", help="MD timestep in picoseconds", default=0.002, type=float)
    parser.add_argument("-ff", "--forcefield", help="Force field to use (can be charmm36, amber14, amoeba, or amber19 on newer openmm versions", default="charmm36")
    parser.add_argument("-i", "--integrator", help="Numerical integration scheme. Can be one of: VariableLangevinIntegrator, LangevinMiddleIntegrator", default="LangevinMiddleIntegrator")
    parser.add_argument("-f", "--friction", help="Friction coefficient which couples the system to the heat bath. Units 1/picosecond.", default=1, type=float)
    parser.add_argument("-t", "--temperature", help="Simulation temperature in Kelvin.", default=300, type=float)
    parser.add_argument("-nomin", "--no_minimise", help="Skip minimisation", action="store_true")
    parser.add_argument("-notest", "--no_testrun", help="Skip test run", action="store_true")
    parser.add_argument("-fb", "--fix_backbone", help="Fix the backbone", action="store_true")
    parser.add_argument("-c", "--constraints", help="Constraints to apply. Possible values: None, HBonds. Constraining HBonds can make simulations run faster, at the cost of some accuracy. You can't constrain HBonds if you're also fixing the backbone, for obvious reasons.", default="Default")
    parser.add_argument("-solv", "--solvent", help="Solvent to use. Possible values: tip3p, tip4pew, spce. If no solvent is specified, none will be used.", default=None)
    parser.add_argument("-ns", "--strip_solvent", help="Don't write solvent molecules to the minimised output file", action="store_true")
    parser.add_argument("-ion", "--ionic_strength", help="Ionic strength of the solvent in molar.", default=0.0, type=float)
    parser.add_argument("-p", "--pressure", help="Pressure in bar, coupling via monte carlo barostat If pressure is not specified, no pressure coupling will be used.", default=None)
    parser.add_argument("-step", "--step", help="In production MD, write to output files (traj, thermo) every x steps.", default=1000, type=int)
    parser.add_argument("-th", "--thermo", help="Thermo information output file.", default="thermo.txt")
    parser.add_argument("-nb", "--nonbonded", help="Non-bonded interaction method. Can be: PME, CutoffPeriodic, CutoffNonPeriodic.", default="Default")
    parser.add_argument("-ckpt", "--checkpoint", help="Checkpoint output filename", default="checkpoint.dat")
    parser.add_argument("-q", "--quiet", help="Don't print info to the stdout", action="store_true")
    parser.add_argument("-parm", "--write_params", help="Write simulation params to a (json) file", default="params.json")
          
    args = parser.parse_args()
    minimise = not args.no_minimise
    testrun = not args.no_testrun
    verbose = not args.quiet
    run(args.pdb,
            minimised_structure_out = args.min_out,
            traj_out = args.traj_out,
            max_minimise_iterations = args.minimise_iterations,
            minimise_error=args.minimise_err,
            test_sim_steps=args.test_steps,
            md_steps=args.md_steps,
            md_timestep=args.timestep*picoseconds,
            forcefield=args.forcefield,
            integrator=args.integrator,
            friction_coeff=args.friction/picosecond,
            temperature=args.temperature*kelvin,
            minimise=minimise,
            test_run=testrun,
            fix_backbone=args.fix_backbone,
            constraints=args.constraints,
            solvent=args.solvent,
            strip_solvent=args.strip_solvent,
            ionic_strength=args.ionic_strength*molar,
            pressure=args.pressure*bar,
            step=args.step,
            thermo_out_file=args.thermo,
            non_bonded_method = args.nonbonded,
            checkpoint_output = args.checkpoint,
            verbose = verbose,
            write_params = args.write_params
        )

if __name__ == "__main__":
    entry_point()


# depreceated
# def variable_minimise(pdb, out, max_iterations=1000, errortol=0.001, steps=50):
#     """
#     Minimise structure with openmm, with a variable langevin integrator.
#     Args:
#         pdb: input structure file (pdb or mmcif), a string
#         out: path to output file, a string
#         max_iterations: max iterations of the miminisation to perform, an int
#         errortol: arbitrary openmm value
#     """
#     # platform, prop = get_prop()
#     if ".cif" in pdb or ".mmcif" in pdb:
#         structure = PDBxFile(pdb)
#         writer = PDBxFile
#     else:
#         structure = PDBFile(pdb)
#         writer = PDBFile
#     forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')
#     system = forcefield.createSystem(structure.topology,
#                                      nonbondedMethod=app.NoCutoff,
#                                      nonbondedCutoff=1*nanometer,
#                                      constraints=HBonds)
#     integrator = VariableLangevinIntegrator(300*kelvin, 1/picosecond, errortol)
#     simulation = Simulation(structure.topology, system, integrator)
#     simulation.context.setPositions(structure.positions)
#     simulation.reporters.append(PDBReporter('output.pdb', 5000))
#     simulation.minimizeEnergy(maxIterations=max_iterations)
#     simulation.step(steps)
#     simulation.minimizeEnergy(maxIterations=50)
#     print("Fixed explosion.")
#     writer.writeModel(structure.topology, simulation.context.getState(
#         getPositions=True).getPositions(asNumpy=True),
#         file=open(out, "w"), keepIds=True)


# def test_sim(pdb, minimise=10, steps=50, timestep=0.002*picoseconds):
#     """
#     Run a short test simulation with openmm.
#     Args:
#         pdb: input structure file (pdb or mmcif), a string
#         minimise: number of minimisation seps, an int
#         steps: number of MD steps to run, an int
#         timestep: timestep in openmm units
#     Returns:
#         nothing. The only thing it might do is throw an error!
#     """
#     # platform, prop = get_prop()
#     if ".cif" in pdb or ".mmcif" in pdb:
#         structure = PDBxFile(pdb)
#     else:
#         structure = PDBFile(pdb)
#     forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')
#     system = forcefield.createSystem(structure.topology,
#                                      nonbondedMethod=app.NoCutoff,
#                                      nonbondedCutoff=1*nanometer,
#                                      constraints=HBonds)
#     integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, timestep)
#     simulation = Simulation(structure.topology, system, integrator)
#     simulation.context.setPositions(structure.positions)
#     simulation.minimizeEnergy(maxIterations=minimise)
#     try:
#         simulation.step(steps)
#         print("Test simulation ran successfully.")
#     except OpenMMException:
#         print("Simulation blew up, trying adaptive minimisation...")
#         variable_minimise(pdb, pdb, max_iterations=1000, errortol=0.001)