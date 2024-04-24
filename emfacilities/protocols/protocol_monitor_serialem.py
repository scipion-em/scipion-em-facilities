# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:    
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import os.path
import os

from pathlib import Path
import numpy as np
import json
import pandas as pd
import time

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import STATUS_FINISHED
from pyworkflow import VERSION_1_1

from .protocol_monitor import ProtMonitor, Monitor
from .protocol_monitor_summary import ProtMonitorSummary
from pyworkflow.protocol.params import PointerParam



class ProtMonitorSerialEm(ProtMonitor):
    """ Stores values for SerialEM to read ,
        out of the ranges defined
    """
    _label = 'monitor SerialEm'
    _lastUpdateVersion = VERSION_1_1

    def __init__(self, **kwargs):
        ProtMonitor.__init__(self, **kwargs)
        self.fileDir = ''
        self.filePath = ''
        self.data = pd.DataFrame()

    def _defineParams(self, form):  

        form.addSection('SerialEm')
        form.addParam('monitorProt', PointerParam,
                      pointerClass=ProtMonitorSummary,
                      label='monitor summary data protocol',
                      help="connnect with the monitor summary to"
                           "extract the necesary info")
        form.addParam('filesPath', params.PathParam,
                      label="SerialEM File directory",
                      help="Directory to store the SerialEM txt.\n\n"
                           "The path can also contain wildcards to select"
                           "from several folders. \n\n"
                           "Examples:\n"
                           "  ~/project/data/day??_files/\n"
                           "Each '?' represents one unknown character\n\n"
                           "  ~/project/data/day*_files/\n"
                           "'*' represents any number of unknown characters\n\n"
                           "  ~/project/data/day##_files/\n"
                           "'##' represents two digits that will be used as "
                           "file ID\n\n"
                           "NOTE: wildcard characters ('*', '?', '#') "
                           "cannot appear in the actual path.)")
        
        form.addParam('maxGlobalShift', params.FloatParam, default=1000,
                      label="maxGlobalShift in a micrograph allowed",
                      help="maximun Global Shift between frames of a "
                            "movie allowed, if parameter is surpassed"
                            "writes 1 to the file")
        form.addParam('maxFrameShift', params.FloatParam, default=100,
                      label="Max Frame Shift allowed ",
                      help="maximun Frame Shift between two consecutive"
                            "Frames of a movie, if parameter is surpassed"
                            "writes 1 to the file")
        form.addParam('maxDefocusU', params.FloatParam, default=0.0,
                      label="Max Defocus U allowed",
                      help="maximun Defocus U (Amstrong) allowed in a micrograph"
                            "if the parameter is surpassed writes 1 "
                            "to the file")
        form.addParam('maxDefocusV', params.FloatParam, default=0.0,
                      label="Max Defocus V allowed",
                      help="maximun Defocus V (Amstrong) allowed in a micrograph"
                            "if the parameter is surpassed writes 1"
                             "to the file")
        
        form.addParam('thresholdShift', params.FloatParam, default=10,
                      label="threshold shift allowed",
                      help="allow to surpass the treshold stablished")
        form.addParam('thresholdDefocus', params.FloatParam, default=5,
                      label="threshold defocus allowed",
                      help="allow to surpass the treshold stablished")
    

    # --------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    # --------------------------- STEPS functions ----------------------------
    def monitorStep(self):
        
        self.fileDir = self.monitorProt.get()._getExtraPath()
        self.filePath = Path(str(self.filesPath)) / "serialEM.csv"
        
        data_dict = {'maxGlobalShift': 0, 'maxFrameShift': 0, 'maxDefocusU': 0, 'maxDefocusV': 0}
        self.data = pd.DataFrame(data_dict, index=[0])

        def readFile():
            original_path = Path(self.fileDir)

            # Get the file paths
            file_mic = original_path / "tmp" / "defocus.txt"
            file_phase = original_path / "tmp" / "phase.txt"

            try:
                with open(file_mic, 'r') as mic_file:
                    defocus_data = json.load(mic_file)

                # Load phases from file_phase
                with open(file_phase, 'r') as phase_file: 
                    phase_data = json.load(phase_file)
            except:
                print("waiting for files to be created")
                return None,None

            return phase_data, defocus_data
        
        def checkFile():

            all_phases, defocus_values = readFile()
            while defocus_values== None:
                all_phases, defocus_values = readFile()
                time.sleep(60)
            for mic,values in defocus_values.items():
                for defocus_U, defocus_V in values:
                    print(f"Defocus_U: {defocus_U}, Defocus_V: {defocus_V}")
                    if defocus_U >= self.maxDefocusU:
                        self.data['maxDefocusU'] = 1

                    if defocus_V >= self.maxDefocusV:
                        self.data['maxDefocusV'] = 1

            threshold=0
            for mic,phase_list in all_phases.items():
                # Define arrays to store shift values for X and Y
                shiftArrayX = phase_list[0]  # X shifts are at the first position
                shiftArrayY = phase_list[1] 

                # Compute frame shifts for X and Y
                frameShiftX = np.abs(np.diff(shiftArrayX))
                frameShiftY = np.abs(np.diff(shiftArrayY))

                maxShiftM = max(shiftArrayX[0]-shiftArrayX[-1], shiftArrayY[0]-shiftArrayY[-1])

                maxShiftBetweenFrames = max(np.max(frameShiftX), np.max(frameShiftY))

                if maxShiftM >= self.maxGlobalShift:
                    if threshold > self.thresholdshift:
                        self.data['maxGlobalShift'] = 1

                if maxShiftM < 0:
                    self.data['maxGlobalShift'] = -1

                if maxShiftBetweenFrames >= self.maxFrameShift:
                    self.data['maxFrameShift'] = 1


            self.data.to_csv(self.filePath, sep='\t', index=False)
        
        while self.monitorProt.get().getStatus() != STATUS_FINISHED:
            checkFile()
            time.sleep(60)

        return None
        


