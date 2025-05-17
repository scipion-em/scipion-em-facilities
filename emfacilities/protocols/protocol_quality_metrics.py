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
import time
import re
from datetime import datetime
import pandas as pd
import numpy as np

import pyworkflow.protocol.params as params
from pyworkflow.protocol import getUpdatedProtocol
from pyworkflow import UPDATED, NEW, VERSION_3_0
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pyworkflow.protocol.constants import STATUS_NEW

from pwem.protocols import ProtImportImages
from pwem.protocols import EMProtocol
from pwem import Domain
import pwem.objects as emobj

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
        self.lastObjCreationTime = {}
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
        self.sleepCheck()

    def sleepCheck(self):
        if not self.finished:
            sleepTime = self.getTimeOutInSeconds(self.samplingInterval.get())
            self.info('Sleeping interval for %d seconds' % sleepTime)
            time.sleep(sleepTime)

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
        self.info('Protocol is active status: %s' % status)
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
            if protName in self.lastObjCreationTime:
                where = 'creation>"' + str(self.lastObjCreationTime[protName]) + '"'

            for cls, action in quality_actions.items():
                if isinstance(prot, cls):
                    self.info("Launching action for protocol %s where %s" % (protName, where))
                    result, df, last_creation_time = action(prot, where)
                    outputs_dicts[protName] = result
                    if last_creation_time: # If no new entries are added the lastCreationTime is not changed
                        self.lastObjCreationTime[protName] = last_creation_time

                    if self.dataFrame.empty:
                        self.info('Creating the first dataFrame')
                        self.dataFrame = df.copy()
                    else:
                        if not df.empty:
                            self.info('Updating entries row-wise')
                            # print(df)
                            # Ensure all columns exist in main DataFrame
                            for col in df.columns:
                                if col != 'movieId' and col not in self.dataFrame.columns:
                                    self.dataFrame[col] = pd.NA
                            # Set movieId as index for efficient lookup
                            self.dataFrame.set_index('movieId', inplace=True)
                            df.set_index('movieId', inplace=True)

                            for movieId, row in df.iterrows():
                                if movieId in self.dataFrame.index:
                                    for col in df.columns:
                                        if pd.notna(row[col]):
                                            self.dataFrame.at[movieId, col] = row[col]
                                else:
                                    # Add new row if movieId not in main DataFrame
                                    self.dataFrame.loc[movieId] = row

                            # Reset index if needed
                            self.dataFrame.reset_index(inplace=True)
                        else:
                            self.info('dataset empty')
                    break

        # Save the current dataframe to CSV after processing
        csv_path = os.path.join(self._getExtraPath(), 'monitor_data.csv')
        self.dataFrame.to_csv(csv_path, index=False)
        self.info(self.dataFrame)
        dict_str = '\n'.join(f'{k}: {v}' for k, v in outputs_dicts.items())
        self.summaryVar.set(dict_str)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        return []

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
    lastCreationTimes = []

    for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
        outSet.load()
        outSet.loadAllProperties()
        pixel_size = round(outSet.getSamplingRate(), 2)
        size = outSet.getSize()
        outputs[outName] = size

        for item in outSet.iterItems(orderBy='creation', direction='ASC', where=where):
            objId = item.getObjId()
            objCreationTime = item.getObjCreation()
            entry = {
                'movieId': objId,
                'movieName': item.getBaseName(),
                'creationTime': objCreationTime,
                'magnification': magnification,
                'pixelSize': pixel_size,  # in Å/pixel
                'voltage': voltage,
                'sphericalAberration': spherical_ab,
                'dosePerFrame': dose_per_frame
            }
            entries.append(entry)

        outSet.close()

        if entries:
            lastCreationTimes.append(objCreationTime)

    df = pd.DataFrame(entries)

    if lastCreationTimes:
        lastCreationTime = lastCreationTimes[-1]
    else:
        lastCreationTime = None

    return outputs, df, lastCreationTime

def extractDoseAnalysis(prot, where):
    percentage_th = prot.percentage_threshold.get()
    outputs = {}
    entries = []
    lastCreationTimes = []

    for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
        outSet.load()
        outSet.loadAllProperties()
        size = outSet.getSize()
        outputs[outName] = size

        boolPassDoseAnalysis = True
        if 'Discarded' in outName:
            boolPassDoseAnalysis = False

        for item in outSet.iterItems(orderBy='creation', direction='ASC', where=where):
            objId = item.getObjId()
            objCreationTime = item.getObjCreation()
            entry = {
                'movieId': objId,
                'movieName': item.getBaseName(),
                'boolPassDoseAnalysis': boolPassDoseAnalysis,
                'thresholdPercentageDiff': percentage_th,
                'diffDosePerAngstrom2': round(item._DIFF_TO_DOSE_PER_ANGSTROM2.get(), 2),
                'meanDosePerAngstrom2': round(item._MEAN_DOSE_PER_ANGSTROM2.get(), 2),
                'stdDosePerAngstrom2': round(item._STD_DOSE_PER_ANGSTROM2.get(), 2)
            }
            entries.append(entry)

        outSet.close()
        if entries:
            lastCreationTimes.append(objCreationTime)

    df = pd.DataFrame(entries)

    if lastCreationTimes:
        lastCreationTime = lastCreationTimes[-1]
    else:
        lastCreationTime = None

    return outputs, df, lastCreationTime

def extractMaxShift(prot, where):
    max_frame_shift_th = prot.maxFrameShift.get()
    max_global_shift_th = prot.maxMovieShift.get()
    outputs = {}
    movieEntries = {}  # indexed by movieId
    lastCreationTimes = []

    for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
        outSet.load()
        outSet.loadAllProperties()
        size = outSet.getSize()
        outputs[outName] = size

        boolPassMaxShift = 'Discarded' not in outName

        if isinstance(outSet, emobj.SetOfMovies):
            pixel_size = round(outSet.getSamplingRate(), 2)
            for item in outSet.iterItems(orderBy='creation', direction='ASC', where=where):
                x_shifts = item._alignment._xshifts
                y_shifts = item._alignment._yshifts
                max_frame_shift, max_movie_shift = calculateMaxShift(x_shifts, y_shifts, pixel_size)
                movieId = item.getObjId()
                objCreationTime = item.getObjCreation()

                entry = {
                    'movieId': movieId,
                    # 'movieName': item.getBaseName(),
                    'boolPassMaxShift': boolPassMaxShift,
                    'thresholdMaxFrameShift': max_frame_shift_th,
                    'thresholdMaxMovieShift': max_global_shift_th,
                    'maxFrameShift': round(max_frame_shift, 2),
                    'maxMovieShift': round(max_movie_shift, 2),
                    'accumMotionTotal': None,
                    'accumMotionEarly': None,
                    'accumMotionLate': None
                }
                movieEntries[movieId] = entry
        else:
            for item in outSet.iterItems(orderBy='creation', direction='ASC', where=where):
                movieId = item.getObjId()
                movieEntries[movieId]['micName'] = item.getBaseName() # MicName
                movieEntries[movieId]['accumMotionTotal'] = round(item._rlnAccumMotionTotal.get(), 2)
                movieEntries[movieId]['accumMotionEarly'] = round(item._rlnAccumMotionEarly.get(), 2)
                movieEntries[movieId]['accumMotionLate'] = round(item._rlnAccumMotionLate.get(), 2)

        outSet.close()
        if movieEntries:
            lastCreationTimes.append(objCreationTime)

    df = pd.DataFrame(movieEntries.values())

    if lastCreationTimes:
        lastCreationTime = lastCreationTimes[-1]
    else:
        lastCreationTime = None

    return outputs, df, lastCreationTime

def calculateMaxShift(x_shifts, y_shifts, pixel_size):
    """
    Calculates max frame-to-frame shift and max movie shift range from absolute shifts.

    Parameters:
        x_shifts (list): Absolute X shifts per frame.
        y_shifts (list): Absolute Y shifts per frame.
        pixel_size (float): Pixel size in Å/pixel.

    Returns:
        max_frame_shift (float): Max per-frame shift (Å).
        max_movie_shift (float): Max range across movie (Å).
    """
    x_shifts = np.asarray(x_shifts)
    y_shifts = np.asarray(y_shifts)
    # Frame-wise shifts (difference between consecutive absolute positions)
    frame_dx = np.diff(x_shifts)
    frame_dy = np.diff(y_shifts)
    frame_shifts = np.sqrt(frame_dx**2 + frame_dy**2)
    max_frame_shift = np.max(frame_shifts) * pixel_size
    # Movie-wise range (max difference in each direction)
    range_x = np.max(x_shifts) - np.min(x_shifts)
    range_y = np.max(y_shifts) - np.min(y_shifts)
    max_movie_shift = max(range_x, range_y) * pixel_size

    return max_frame_shift, max_movie_shift

def extractTiltAnalaysis(prot, where):
    mean_corr_th = round(prot.meanCorr_threshold.get(), 2)
    std_corr_th = round(prot.stdCorr_threshold.get(), 2)

    outputs = {}
    entries = []
    lastCreationTimes = []

    for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
        outSet.load()
        outSet.loadAllProperties()
        size = outSet.getSize()
        outputs[outName] = size

        boolPassTiltAnalysis = True
        if 'Discarded' in outName:
            boolPassTiltAnalysis = False

        for item in outSet.iterItems(orderBy='creation', direction='ASC', where=where):
            objId = item.getObjId()
            objCreationTime = item.getObjCreation()
            entry = {
                'movieId': objId,
                'micName': item.getBaseName(),
                'boolPassTiltAnalysis': boolPassTiltAnalysis,
                'thresholdMeanCorrelation': mean_corr_th,
                'thresholdStdCorrelation': std_corr_th,
                'tiltMeanCorrelation': round(item._tilt_mean_corr.get(), 2),
                'tiltStdCorrelation': round(item._tilt_std_corr.get(), 2),
                'tiltImage': item._tilt_psds_image._filename.get()
            }
            entries.append(entry)

        outSet.close()
        if entries:
            lastCreationTimes.append(objCreationTime)

    df = pd.DataFrame(entries)

    if lastCreationTimes:
        lastCreationTime = lastCreationTimes[-1]
    else:
        lastCreationTime = None

    return outputs, df, lastCreationTime

def extractCTFConsensus(prot, where):
    max_defocus_th = prot.maxDefocus.get()
    min_defocus_th = prot.minDefocus.get()
    astigmatism_percentage_th = round(prot.astigmatismPer.get(),2)
    resolution_th = prot.resolution.get()
    consensus_resolution_th = prot.minConsResol.get() if prot.calculateConsensus.get() else None

    outputs = {}
    entries = []
    lastCreationTimes = []

    for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
        outSet.load()
        outSet.loadAllProperties()
        size = outSet.getSize()
        outputs[outName] = size

        boolPassCTFConsensus = 'Discarded' not in outName

        if isinstance(outSet, emobj.SetOfCTF):
            for item in outSet.iterItems(orderBy='creation', direction='ASC', where=where):
                objId = item.getObjId()
                objCreationTime = item.getObjCreation()
                entry = {
                    'movieId': objId,
                    'boolPassCTFConsensus': boolPassCTFConsensus,
                    'thresholdMaxDefocus': max_defocus_th,
                    'thresholdMinDefocus': min_defocus_th,
                    'thresholdAstigmatismPercentage': astigmatism_percentage_th,
                    'thresholdResolution': resolution_th,
                    'thresholdConsensusResolution': consensus_resolution_th,
                    'defocusU': round(item._defocusU.get(), 2),
                    'defocusV': round(item._defocusV.get(), 2),
                    'defocusRatio': round(item._defocusRatio.get() ,2),
                    'defocusAngle': round(item._defocusAngle.get(), 2),
                    'astigmatismPercentage': round(item._astigmatismPercentage.get(), 2),
                    'resolution': round(item._resolution.get(), 2),
                    'fitQuality': round(item._fitQuality.get(), 2),
                    'IceRingDensity': round(item._rlnIceRingDensity.get(), 2),
                    'consensusResolution': round(item._consensus_resolution.get(), 2) if consensus_resolution_th else np.NaN,
                    'psdFile': item._psdFile.get(),
                }
                entries.append(entry)

            if entries:
                lastCreationTimes.append(objCreationTime)

        outSet.close()

    df = pd.DataFrame(entries)

    if lastCreationTimes:
        lastCreationTime = lastCreationTimes[-1]
    else:
        lastCreationTime = None

    return outputs, df, lastCreationTime

def extractMiffi(prot, where):
    outputs = {}
    entries = []
    lastCreationTimes = []

    for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
        outSet.load()
        outSet.loadAllProperties()
        size = outSet.getSize()
        outputs[outName] = size

        boolPassMiffi = True
        if 'Discarded' in outName:
            boolPassMiffi = False

        for item in outSet.iterItems(orderBy='creation', direction='ASC', where=where):
            objId = item.getObjId()
            objCreationTime = item.getObjCreation()
            entry = {
                'movieId': objId,
                'micName': item.getBaseName(),
                'boolPassMiffi': boolPassMiffi,
                'miffiLabel': item._miffi_label.get()
            }
            entries.append(entry)

        outSet.close()
        if entries:
            lastCreationTimes.append(objCreationTime)

    df = pd.DataFrame(entries)

    if lastCreationTimes:
        lastCreationTime = lastCreationTimes[-1]
    else:
        lastCreationTime = None

    return outputs, df, lastCreationTime

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


