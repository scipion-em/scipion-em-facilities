# **************************************************************************
# *
# * Authors: Daniel Marchan (da.marchan@cnb.csic.es)
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
from datetime import datetime, timedelta
import time
import copy
import re

from pwem.protocols import EMProtocol
from pwem.objects import SetOfImages, Set

from pyworkflow import VERSION_3_0
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow import UPDATED, NEW
from pyworkflow.protocol.constants import STATUS_NEW



OUTPUT = "outputSet"

class ProtDataCounter(EMProtocol):
    """
    Protocol to make a subset of images from the original one.
    Waits until certain number of images is prepared and then send them to output.
    It can be done in 3 ways:
        - If simple mode: once the number of items is reached, a setOfImages is returned and
            the protocol finishes (ending the streaming from this point).
    The protocol will accept Micrographs, Particles and any kind of object that inherits from Image base class.
    """
    _label = 'data counter'
    _devStatus = NEW
    _lastUpdateVersion = VERSION_3_0
    _possibleOutputs = {OUTPUT: SetOfImages}


    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputImages', params.PointerParam, pointerClass='SetOfImages',
                      label="Input images", important=True)
        form.addParam('outputSize', params.IntParam, default=10000,
                      label='Output size',
                      help='How many images need to be on input to '
                           'create output set.')

        form.addSection(label='Stream of data timer')
        form.addParam('delay', params.IntParam, default=10, label="Delay (sec)",
                      validators=[params.GT(2, "must be larger than 3sec.")],
                      help="Delay in seconds before checking new output")
        form.addParam('boolTimer', params.BooleanParam, default=False,
                      label='Use timer?',
                      help='Select YES if you want to use a timer to close the stream of data once the timer is consumed.\n'
                           'If NO is selected, normal functionality.\n'
                           'If the output size is achieved first it will be closed as well.')
        form.addParam('timeout', params.StringParam, default="1h",
                      condition="boolTimer==%d" % True,
                      label='Time to wait:',
                      help='Time in seconds that the protocol will remain '
                           'running. A correct format is an integer number in '
                           'seconds or the following syntax: {days}d {hours}h '
                           '{minutes}m {seconds}s separated by spaces '
                           'e.g: 1d 2h 20m 15s,  10m 3s, 1h, 20s or 25.')


# --------------------------- INSERT steps functions -------------------------
    def _insertAllSteps(self):
        self.initializeParams()
        self._insertFunctionStep(self.createOutputStep,
                                 prerequisites=[], wait=True, needsGPU=False)

    def createOutputStep(self):
        self._closeOutputSet()

    def initializeParams(self):
        self.finished = False
        # Important to have both:
        self.insertedIds = []   # Contains images that have been inserted in a Step (checkNewInput).
        self.processedIds = [] # Ids to be output
        self.isStreamClosed = self.inputImages.get().isStreamClosed()
        # Contains images that have been processed in a Step (checkNewOutput).
        self.inputFn = self.inputImages.get().getFileName()
        self._inputClass = self.inputImages.get().getClass()
        self._inputType = self.inputImages.get().getClassName().split('SetOf')[1]
        self._baseName = '%s.sqlite' % self._inputType.lower()
        self.limitReach = False
        self.timerOut = False
        self.timeoutSecs = self.getTimeOutInSeconds(self.timeout.get())
        self.lastTimeCheckTimer = datetime.now() # Timer

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
        self.lastCheck = getattr(self, 'lastCheck', datetime.now())
        mTime = datetime.fromtimestamp(os.path.getmtime(self.inputFn))
        self.debug('Last check: %s, modification: %s'
                    % (pwutils.prettyTime(self.lastCheck),
                        pwutils.prettyTime(mTime)))
        # If the input.sqlite have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime and self.insertedIds:  # If this is empty it is due to a static "continue" action or it is the first round
            return None

        inputSet = self._loadInputSet(self.inputFn)
        inputSetIds = inputSet.getIdSet()
        newIds = [idImage for idImage in inputSetIds if idImage not in self.insertedIds]

        self.lastCheck = datetime.now()
        self.isStreamClosed = inputSet.isStreamClosed()
        inputSet.close()

        outputStep = self._getFirstJoinStep()

        if self.isContinued() and not self.insertedIds:  # For "Continue" action and the first round
            doneIds, _ = self._getAllDoneIds()
            skipIds = list(set(newIds).intersection(set(doneIds)))
            newIds = list(set(newIds).difference(set(doneIds)))
            self.info("Skipping Images with ID: %s, seems to be done" % skipIds)
            self.insertedIds = doneIds  # During the first round of "Continue" action it has to be filled

        if newIds and not self.limitReach and not self.timerOut:
            fDeps = self._insertNewImageSteps(newIds)
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()

    def _checkNewOutput(self):
        doneListIds, currentOutputSize = self._getAllDoneIds()
        processedIds = copy.deepcopy(self.processedIds)
        newDone = [imageId for imageId in processedIds if imageId not in doneListIds]
        allDone = len(doneListIds) + len(newDone)
        limitOutputSize = self.outputSize.get()
        maxSize = self._loadInputSet(self.inputFn).getSize()
        self.limitReach = allDone >= limitOutputSize

        # We have finished when there is not more input images
        # (stream closed) or when the limit of output size is met
        self.finished = (self.isStreamClosed and allDone == maxSize) or (self.limitReach or self.timerOut)

        if not self.finished and not newDone:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        inputSet = self._loadInputSet(self.inputFn)
        outputSet = self._loadOutputSet(self._inputClass, self._baseName)

        if currentOutputSize < limitOutputSize:
            for imageId in newDone:
                image = inputSet.getItem("id", imageId).clone()
                outputSet.append(image)
                currentOutputSize += 1
                if currentOutputSize == limitOutputSize:
                    self.finished = True
                    break # We have reach the limit for the outputSize

            streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

            self._updateOutputSet(OUTPUT, outputSet, streamMode)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

        self._store()

    def _loadInputSet(self, inputFn):
        self.debug("Loading input db: %s" % inputFn)
        inputSet = self._inputClass(filename=inputFn)
        inputSet.loadAllProperties()
        return inputSet

    def _loadOutputSet(self, SetClass, baseName):
        setFile = self._getPath(baseName)

        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        inputs = self.inputImages.get()
        outputSet.copyInfo(inputs)

        return outputSet

    def _insertNewImageSteps(self, newIds):
        """ Insert steps to register new images (from streaming)
        Params:
            newIds: input images ids to be processed
        """
        deps = []
        stepId = self._insertFunctionStep(self.registerStep, newIds, needsGPU=False,
                                          prerequisites=[])
        deps.append(stepId)
        self.insertedIds.extend(newIds)

        return deps

    def registerStep(self, newIds):
        self.info('Registering the %d new images' %len(newIds))
        for imageId in newIds:
            self.processedIds.append(imageId)

        if not self.isStreamClosed:
            self.delayRegister()
            if self.boolTimer.get():
                self.info('Using timer:')
                self.timerStep()

    def delayRegister(self):
        delay = self.delay.get()
        self.info('Sleeping for delay time: %d' %delay)
        time.sleep(delay)

    def timerStep(self):
        endTime = self.lastTimeCheckTimer + timedelta(seconds=self.timeoutSecs)
        now = datetime.now()
        remainingTime = (endTime - now).total_seconds()

        if remainingTime <= 0:
            self.timerOut = True
            self.info("  timer is consumed terminating protocol.")
        else:
            self.info(f"  remaining time: {int(remainingTime)} seconds.")
            self.timeoutSecs = int(remainingTime)

        # Update the last time check
        self.lastTimeCheckTimer = now

    # ------------------------- UTILS functions --------------------------------
    def _getAllDoneIds(self):
        doneIds = []
        sizeOutput = 0

        if hasattr(self, OUTPUT):
            sizeOutput = self.outputSet.getSize()
            doneIds.extend(list(self.outputSet.getIdSet()))

        return doneIds, sizeOutput

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

    def _summary(self):
        pass