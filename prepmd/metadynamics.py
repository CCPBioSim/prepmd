#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run metadynamics simulations - used in run.py
"""

from openmm.unit import *
import numpy as np
import MDAnalysis.analysis.rms as rms

class StopSimulation(Exception):
    """ Exception raised in order to stop simulation """
    pass

class RMSDIncrease(Exception):
    pass

class RMSDReporter(object):
    def __init__(self, file, reportInterval, ref_positions, threshold_nm = 0.17, consecutive_frames=2):
        self._out = open(file, 'w')
        self._reportInterval = reportInterval
        self._ref_positions = ref_positions.value_in_unit(nanometers)
        self._threshold_nm = threshold_nm
        self.max_consecutive_frames = consecutive_frames
        self.current_consecutive_frames = 0
        self.rmsd_log = []
        self.reports = 0

    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return {'steps': steps, 'periodic': None, 'include':['positions']}

    def report(self, simulation, state):
        positions = state.getPositions(asNumpy=True).value_in_unit(nanometers)
        # i think this measurement doesn't align whereas rmsdforce does
        #rmsd = np.sqrt(((positions - self._ref_positions)**2).sum(-1).mean())
        report_rmsd = rms.rmsd(positions, self._ref_positions, center=True, superposition=True)
        print("rmsd = "+str(report_rmsd), end=" ")
        self.rmsd_log.append(report_rmsd)
        self._out.write(str(report_rmsd)+"\n")
        if report_rmsd < self._threshold_nm:
            self.current_consecutive_frames += 1
        else:
            self.current_consecutive_frames = 0
        if self.current_consecutive_frames > self.max_consecutive_frames:
            print(". RMSD below threshold value, stopping...")
            raise StopSimulation()
        if self.reports%5 == 0 and self.reports > 0:
            if np.mean(np.gradient(self.rmsd_log)) > 0:
                raise RMSDIncrease
        self.reports += 1

def get_representative_rmsd_frames(rmsd_log, num_frames):
    num_frames -= 2
    frames = [0]
    change = 0
    for frame in range(len(rmsd_log)-1):
        change += abs(rmsd_log[frame+1] - rmsd_log[frame])
    change_per_frame = change / num_frames
    change = 0
    for frame in range(len(rmsd_log)-1):
        change += abs(rmsd_log[frame+1] - rmsd_log[frame])
        if change > change_per_frame:
            frames.append(frame)
            change -= change_per_frame
    frames.append(len(rmsd_log)-1)
    return frames