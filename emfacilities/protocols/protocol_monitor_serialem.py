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
import sys

from pathlib import Path
import numpy as np
import json
import pandas as pd
import time

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import STATUS_FINISHED, STATUS_FAILED, STATUS_ABORTED
from pyworkflow import VERSION_1_1

from .protocol_monitor import ProtMonitor, Monitor
from .protocol_monitor_summary import ProtMonitorSummary
from pyworkflow.protocol.params import PointerParam

STATUS_STOP = [STATUS_FINISHED, STATUS_FAILED, STATUS_ABORTED]

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
        self.checkStatus = None

    def _defineParams(self, form):  

        form.addSection('SerialEm')
        form.addParam('monitorProt', PointerParam,
                      pointerClass=ProtMonitorSummary,
                      label='Monitor summary data protocol',
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
                           "'*'self.monitorProt.get().getStatus() represents any number of unknown characters\n\n"
                           "  ~/project/data/day##_files/\n"
                           "'##' represents two digits that will be used as "
                           "file ID\n\n"
                           "NOTE: wildcard characters ('*', '?', '#') "
                           "cannot appear in the actual path.)")
        
        form.addParam('maxGlobalShift', params.FloatParam, default=20,
                      label="MaxGlobalShift in a micrograph allowed",
                      help="maximun Global Shift between frames of a "
                            "movie allowed, if parameter is surpassed"
                            "writes 1 to the file")
        form.addParam('maxFrameShift', params.FloatParam, default=5,
                      label="Max Frame Shift allowed ",
                      help="maximun Frame Shift between two consecutive"
                            "Frames of a movie, if parameter is surpassed"
                            "writes 1 to the file")
        form.addParam('maxDefocus', params.FloatParam, default=35000,
                      label="Max Defocus allowed",
                      help="insert positive values maximun Defocus U (Angstrong) allowed in a micrograph"
                            "if the parameter is surpassed writes 1 "
                            "to the file")
        form.addParam('minDefocus', params.FloatParam, default=5000,
                      label="Min Defocus allowed",
                      help="insert positive values maximun DUpdateefocus V (Angstrong) allowed in a micrograph"
                            "if the parameter is surpassed writes -1"
                             "to the file")
        form.addParam('astigmatism', params.FloatParam, default=1.05,
                      label="Astigmatism tolerance",
                      help="allow to surpass the treshold stablished")
        form.addParam('activatewriting', params.BooleanParam, default=True,
                      label="Activate writing in file",
                      help="writes in the file inside the folder")
        form.addParam('afis', params.FloatParam, default=5,
                      label="AFIS",
                      help="number of micrographs per photo taken")
        form.addParam('setafis', params.FloatParam, default=2,
                      label="Number of AFIS",
                      help="Number of AFIS consecutively surpassed to activate"
                            "astigmatism correction")
    # --------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    # --------------------------- STEPS functions ----------------------------
    def monitorStep(self):
        
        self.fileDir = self.monitorProt.get()._getExtraPath()
        self.filePath = Path(str(self.filesPath)) / "serialEM.csv"
        
        data_dict = {'maxGlobalShift': 0, 'maxFrameShift': 0, 'rangeDefocus': 0, 'astigmatism': 0}
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
            afiscount=0 
            micount=0 # count the number of mics
            evaluated_afis = set()
            for mic,values in defocus_values.items():
                for defocus_U, defocus_V in values:
                    afis_key = f"{mic}_{defocus_U}_{defocus_V}"  # Create a unique key for the AFIS values
                    
                    if defocus_U >= self.maxDefocus.get() and  afis_key not in evaluated_afis: 
                        self.data['maxDefocusU'] = 1
                        print("Defocus_U exceeded range",defocus_U)
                        
                    if defocus_U <= self.minDefocus.get() and  afis_key not in evaluated_afis:
                        self.data['maxDefocusU'] = -1
                        print("Defocus_U exceeded range",defocus_U)
                        
                    if defocus_U >= defocus_V and  afis_key not in evaluated_afis:                      
                        ratio = defocus_U / defocus_V
                    else:
                        ratio = defocus_V / defocus_U
                          
                    if ratio >= self.astigmatism.get() and  afis_key not in evaluated_afis :
                        afiscount=afiscount+1
                        if self.afis.get()*self.setafis.get()*0.2 >= afiscount :
                            self.data['astigmatism'] = 1
                            print("astigmatism exceeded range",ratio)
    
                    if afis_key not in evaluated_afis:  # Check if the AFIS values have not been evaluated before
                        evaluated_afis.add(afis_key) 
                    
                    micount=micount+1

                    if micount> self.afis.get()*self.setafis.get():
                        afiscount=0
                        micount=0

            for mic,phase_list in all_phases.items():
                # Define arrays to store shift values for X and Y
                shiftArrayX = phase_list[0]  # X shifts are at the first position
                shiftArrayY = phase_list[1] 

                # Compute frame shifts for X and Y
                frameShiftX = np.abs(np.diff(shiftArrayX))
                frameShiftY = np.abs(np.diff(shiftArrayY))

                maxShiftM = max(shiftArrayX[0]-shiftArrayX[-1], shiftArrayY[0]-shiftArrayY[-1])

                maxShiftBetweenFrames = max(np.max(frameShiftX), np.max(frameShiftY))
                
                if maxShiftM >= self.maxGlobalShift.get() :
                    if self.activatewriting.get():
                        self.data['maxGlobalShift'] = 1
                        print("maxGlobalShift exceeded range",maxShiftM)

                if maxShiftM < -self.maxGlobalShift.get()  :
                    self.data['maxGlobalShift'] = -1
                    print("maxGlobalShift exceeded range",maxShiftM)

                if maxShiftBetweenFrames >= self.maxFrameShift.get() :
                    self.data['maxFrameShift'] = 1
                    print("maxFrameShift exceeded range",maxShiftBetweenFrames)
                

            self.data.to_csv(self.filePath, sep='\t', index=False)
        
        while self.checkStatus not in STATUS_STOP:
            checkFile()
            self.checkStatus = self.monitorProt.get().getStatus()
            self.info('Checking protocol status: %s' % self.checkStatus)

            if self.checkStatus in STATUS_STOP:
                self.info('Finishing protocol')
                continue

            time.sleep(60)
            sys.stdout.flush()

        return None
        


