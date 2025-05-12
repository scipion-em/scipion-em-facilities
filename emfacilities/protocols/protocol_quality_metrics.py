# **************************************************************************
# *
# * Authors:     Daniel Marchan Torres (da.marchan@cnb.csic.es)
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
import os
import sys
import time
import re
from datetime import datetime
import pandas as pd

import pyworkflow.protocol.params as params
from pyworkflow.protocol import getUpdatedProtocol
from pyworkflow import UPDATED, NEW, VERSION_3_0
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pyworkflow.protocol.constants import STATUS_NEW

from pwem.protocols import ProtImportImages
from pwem.protocols import EMProtocol
from pwem import Domain

XMIPP_QUALITY_PROTOCOLS = ['XmippProtMovieDoseAnalysis', 'XmippProtMovieMaxShift', 'XmippProtTiltAnalysis', 'XmippProtCTFConsensus']
MIFFI_QUALITY_PROTOCOLS = ['MiffiProtMicrographs']



class ProtQualityMetrics(EMProtocol):
    """
    """

    _label = "quality metrics"
    _devStatus = NEW
    _lastUpdateVersion = VERSION_3_0


    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputProtocols', params.MultiPointerParam,
                      label="Input protocols", important=True,
                      pointerClass='EMProtocol',
                      help="this protocol/s will be monitorized")

        form.addParam('samplingInterval', params.StringParam, default="60s",
                      label="Sampling Interval:",
                      help='Time in seconds that the protocol will remain '
                           'running. A correct format is an integer number in '
                           'seconds or the following syntax: {days}d {hours}h '
                           '{minutes}m {seconds}s separated by spaces '
                           'e.g: 1d 2h 20m 15s,  10m 3s, 1h, 20s or 25.')

        form.addParam('monitorTime', params.FloatParam, default=34560,
                      label="Total Logging time (min)",
                      help="Log during this interval. 21 days by default")

    # -------------------------- INSERT steps functions -----------------------

    def _insertAllSteps(self):
        self.initializeParams()
        self._insertFunctionStep(self.createOutputStep,
                                 prerequisites=[], wait=True, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------

    def createOutputStep(self):
        self._closeOutputSet()

    def initializeParams(self):
        self.finished = False
        self.lastCreationTime = {}
        self.protsStepDb = [p.getStepsFile() for p in self.getInputProtocols()]
        self.checkingRounds = []
        self.dataFrame = pd.DataFrame()

    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all ctfs
        # to have completed, this can be overriden in subclasses
        # (e.g., in Xmipp 'sortPSDStep')
        return 'createOutputStep'

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def _stepsCheck(self):
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        # Check if there are new images to process from the input set
        if self.finished:
            return

        self.lastCheck = getattr(self, 'lastCheck', datetime.now())
        mTime = self.getModificationTime()

        self.debug('Last check: %s, modification: %s'
                   % (pwutils.prettyTime(self.lastCheck),
                      pwutils.prettyTime(mTime)))

        # If the input.sqlite have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime and self.checkingRounds:  # If this is empty it is due to a static "continue" action or it is the first round
            return None

        self.lastCheck = datetime.now()
        outputStep = self._getFirstJoinStep()

        fDep = self._insertFunctionStep(self.monitorStep, needsGPU=False, prerequisites=[])
        self.checkingRounds.append(fDep)

        if outputStep is not None:
            outputStep.addPrerequisites(fDep)

        self.updateSteps()

    def _checkNewOutput(self):
        prots = [getUpdatedProtocol(p) for p in self.getInputProtocols()]
        status = [p.isActive() for p in prots]
        print(status)

        # We have finished when all protcols have finished
        # (stream closed) or when the limit of output size is met
        self.finished = not any(status)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

        self._store()

    def getModificationTime(self):
        times = []
        for fileDb in self.protsStepDb:
            if os.path.exists(fileDb):
                mTime = datetime.fromtimestamp(os.path.getmtime(fileDb))
                times.append(mTime)

        return max(times)

    def monitorStep(self):
        outputs_dicts = {}
        prots = [getUpdatedProtocol(p) for p in self.getInputProtocols()]
        quality_actions = getImportProtocolAction()
        quality_actions.update(getXmippQualityProtocolsActions())
        quality_actions.update(getMiffiQualityProtocolsActions())
        where = None

        for prot in prots:
            protName = '%s (id=%s)' % (prot.getRunName(), prot.strId())

            if protName in self.lastCreationTime:
                where = 'creation>"' + str(self.lastCreationTime[protName]) + '"'

            for cls, action in quality_actions.items():
                if isinstance(prot, cls):
                    result, df, last_time = action(prot, where)
                    outputs_dicts[protName] = result
                    self.lastCreationTime[protName] = last_time
                    if self.dataFrame.empty:
                        self.info('Creation of the first dataFrame')
                        self.dataFrame = df
                    else:
                        self.info('Dataframe Merge')
                        self.dataFrame = pd.merge(self.dataFrame, df, on='movieId', how='left')
                    break

        print(self.dataFrame)
        dict_str = '\n'.join(f'{k}: {v}' for k, v in outputs_dicts.items())
        self.summaryVar.set(dict_str)
        sleepTime = self.getTimeOutInSeconds(self.samplingInterval.get())
        time.sleep(sleepTime)
        self.info('Sleeping for %d seconds' % sleepTime)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        return []  # no errors

    def _summary(self):
        summary = []
        summary.append(self.summaryVar.get())
        return summary

    def _methods(self):
        return []

    @classmethod
    def worksInStreaming(cls):
        # A monitor protocol always work in streaming
        return True

    # ------------------------------------------ Utils ---------------------------------------------------

    def getInputProtocols(self):
        protocols = []
        for protPointer in self.inputProtocols:
            prot = protPointer.get()
            prot.setProject(self.getProject())
            protocols.append(prot)
        return protocols

    def getTimeOutInSeconds(self, timeOut):
        timeOutFormatRegexList = {r'\d+s': 1, r'\d+m': 60, r'\d+h': 3600,
                                  r'\d+d': 72000}
        try:
            return int(timeOut)
        except Exception:
            seconds = 0
        for regex, secondsUnit in timeOutFormatRegexList.items():
            matchingTimes = re.findall(regex, timeOut)
            for matchTime in matchingTimes:
                seconds += int(matchTime[:-1]) * secondsUnit

        return seconds

def extractImportMovies(prot, where):
    # Common information
    voltage = prot.voltage.get()
    spherical_ab = prot.sphericalAberration.get()
    magnification = prot.magnification.get()
    dose_per_frame = prot.dosePerFrame.get() if prot.dosePerFrame.hasValue() else None

    outputs = {}
    entries = []

    for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
        outSet.load()
        outSet.loadAllProperties()
        pixel_size = round(outSet.getSamplingRate(), 2)
        size = outSet.getSize()
        outputs[outName] = size

        for item in outSet.iterItems(orderBy='creation', direction='ASC', where=where):
            entry = {
                'movieId': item.getObjId(),
                'movieName': item.getBaseName(),
                'magnification': magnification,
                'pixelSize': pixel_size,  # in Å/pixel
                'voltage': voltage,
                'sphericalAberration': spherical_ab,
                'dosePerFrame': dose_per_frame,
                'boolDoseAnalysis': None,
                'thresholdPercentageDiff': None,
                'diffDosePerAngstrom2': None,
                'meanDosePerAngstrom2': None,
                'stdDosePerAngstrom2':None
            }
            entries.append(entry)
            lastCreationTime = item.getObjCreation()

        outSet.close()

    df = pd.DataFrame(entries)

    return outputs, df, lastCreationTime

def extractDoseAnalysis(prot, where):
    percentage_th = prot.percentage_threshold.get()
    outputs = {}
    entries = []
    lastCreationsTimes = []

    for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
        outSet.load()
        outSet.loadAllProperties()
        size = outSet.getSize()
        outputs[outName] = size

        boolPassDoseAnalysis = True
        if 'Discarded' in outName:
            boolPassDoseAnalysis = False

        for item in outSet.iterItems(orderBy='creation', direction='ASC', where=where):
            entry = {
                'movieId': item.getObjId(),
                'movieName': item.getBaseName(),
                'magnification': None,
                'pixelSize': None,  # in Å/pixel
                'voltage': None,
                'sphericalAberration': None,
                'dosePerFrame': None,
                'boolPassDoseAnalysis': boolPassDoseAnalysis,
                'thresholdPercentageDiff': percentage_th,
                'diffDosePerAngstrom2': item._DIFF_TO_DOSE_PER_ANGSTROM2,
                'meanDosePerAngstrom2': item._MEAN_DOSE_PER_ANGSTROM2,
                'stdDosePerAngstrom2': item._STD_DOSE_PER_ANGSTROM2
            }
            entries.append(entry)
            lastCreationTime = item.getObjCreation()

        outSet.close()
        lastCreationsTimes.append(lastCreationTime)

    df = pd.DataFrame(entries)

    return outputs, df, max(lastCreationTime)


def extractMaxShift(prot):
    outputs = {}
    for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
        outSet.load()
        outSet.loadAllProperties()
        # outSetId needs to be compound id to avoid duplicate ids
        outSetId = '%s.%s' % (outSet.getObjId(), prot.getObjId())
        # addObj(outSetId, '', outName, outSet.getSize(), pobj)
        size = outSet.getSize()
        outputs[outName] = size
        outSet.close()

    return outputs

def extractTiltAnalaysis(prot):
    outputs = {}
    for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
        outSet.load()
        outSet.loadAllProperties()
        # outSetId needs to be compound id to avoid duplicate ids
        outSetId = '%s.%s' % (outSet.getObjId(), prot.getObjId())
        # addObj(outSetId, '', outName, outSet.getSize(), pobj)
        size = outSet.getSize()
        outputs[outName] = size
        outSet.close()

    return outputs

def extractCTFConsensus(prot):
    outputs = {}
    for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
        outSet.load()
        outSet.loadAllProperties()
        # outSetId needs to be compound id to avoid duplicate ids
        outSetId = '%s.%s' % (outSet.getObjId(), prot.getObjId())
        # addObj(outSetId, '', outName, outSet.getSize(), pobj)
        size = outSet.getSize()
        outputs[outName] = size
        outSet.close()

    return outputs

def extractMiffi(prot):
    outputs = {}
    for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
        outSet.load()
        outSet.loadAllProperties()
        # outSetId needs to be compound id to avoid duplicate ids
        outSetId = '%s.%s' % (outSet.getObjId(), prot.getObjId())
        # addObj(outSetId, '', outName, outSet.getSize(), pobj)
        size = outSet.getSize()
        outputs[outName] = size
        outSet.close()

    return outputs

def getImportProtocolAction():
    prot_actions = {}
    prot_actions[ProtImportImages] = extractImportMovies

    return prot_actions

def getXmippQualityProtocolsActions():
    actions = [extractDoseAnalysis, extractMaxShift, extractTiltAnalaysis, extractCTFConsensus]
    prot_actions = {}
    for index, classProt in enumerate(XMIPP_QUALITY_PROTOCOLS):
        prot = Domain.importFromPlugin('xmipp3.protocols', classProt)
        if prot is not None:
            prot_actions[prot] = actions[index]
        else:
            print('Problems importing protocol: %s' % classProt)

    return prot_actions

def getMiffiQualityProtocolsActions():
    actions = [extractMiffi]
    prot_actions = {}
    for index, classProt in enumerate(MIFFI_QUALITY_PROTOCOLS):
        prot = Domain.importFromPlugin('miffi.protocols', classProt)
        if prot is not None:
            prot_actions[prot] = actions[index]
        else:
            print('Problems importing protocol: %s' % classProt)

    return prot_actions


