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
import hashlib
import json
import time
import copy
import random

from pyworkflow import VERSION_3_0
from pwem.objects import SetOfImages, Set
import pyworkflow.protocol.params as params

from pwem.protocols import EMProtocol
from pyworkflow import UPDATED, NEW
from pyworkflow.protocol.constants import STATUS_NEW



OUTPUT = "outputSet"

class ProtDataSampler(EMProtocol):
    """
    Protocol to make a subset of images from the original one.
    Waits until certain batch of images is prepared, then it samples a percentage of it and send them to output.
    The protocol will accept Micrographs, Particles, ..., and any kind of object that inherits from Image base class.
    """
    _label = 'data sampler'
    _devStatus = NEW
    _lastUpdateVersion = VERSION_3_0
    _possibleOutputs = {OUTPUT: SetOfImages}


    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputImages', params.PointerParam, pointerClass='SetOfImages',
                      label="Input images", important=True)
        form.addParam('batchSize', params.IntParam, default=10000,
                      label='Batch size',
                      help='How many images need to be on input to '
                           'make the random sampling.')
        form.addParam('samplingProportion', params.FloatParam, default=0.25,
                      label='Sampling proportion',
                      help='What proportion of images need to be output from '
                           'the random sampling (1 means all images and 0 means none).')

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
        self.insertedIds = set()   # Contains images that have been inserted in a Step (checkNewInput).
        self.processedIds = set() # Ids to be register to output
        self.sampleIds = set() # Ids to be output
        self.isStreamClosed = self.inputImages.get().isStreamClosed()
        # Contains images that have been processed in a Step (checkNewOutput).
        self.inputFn = self.inputImages.get().getFileName()
        self._inputClass = self.inputImages.get().getClass()
        self._inputType = self.inputImages.get().getClassName().split('SetOf')[1]
        self._baseName = '%s.sqlite' % self._inputType.lower()

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
        # Always inspect the logical input set. File mtimes are not a valid
        # change detector when the set is backed by PostgreSQL.
        inputSet = self._loadInputSet(self.inputFn)
        inputSetIds = inputSet.getIdSet()

        self.isStreamClosed = inputSet.isStreamClosed()
        inputSet.close()

        outputStep = self._getFirstJoinStep()

        if self.isContinued() and not self.insertedIds:
            doneIds, _ = self._getAllDoneIds()
            self._restoreRuntimeStateFromFinishedSteps(doneIds)
            skipIds = list(set(inputSetIds).intersection(self.insertedIds))
            self.info("Skipping Images with ID: %s, seems to be done" % skipIds)

        newIds = [
            imageId
            for imageId in inputSetIds
            if imageId not in self.insertedIds
        ]

        # Now handle the steps depending on the streaming batch size
        batchSize = self.batchSize.get()
        if len(newIds) < batchSize and not self.isStreamClosed:
            return

        if newIds:
            fDeps = self._insertNewImageSteps(newIds, batchSize)
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()


    def _checkNewOutput(self):
        doneListIds, currentOutputSize = self._getAllDoneIds()
        doneIdSet = set(doneListIds)
        newDone = list(self.sampleIds - doneIdSet)

        inputSet = self._loadInputSet(self.inputFn)
        try:
            inputSetIds = set(inputSet.getIdSet())
            self.finished = (
                self.isStreamClosed
                and inputSetIds.issubset(self.processedIds)
            )

            streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

            if not self.finished and not newDone:
                return

            outputSet = self._loadOutputSet(self._inputClass, self._baseName)

            for imageId in newDone:
                image = inputSet.getItem("id", imageId).clone()
                outputSet.append(image)

            self._updateOutputSet(OUTPUT, outputSet, streamMode)
        finally:
            inputSet.close()

        if self.finished:
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

    def _insertNewImageSteps(self, newIds, batchSize):
        """ Insert steps to register new images (from streaming)
        Params:
            newIds: input images ids to be processed
        """
        deps = []
        # Loop through the image IDs in batches
        for i in range(0, len(newIds), batchSize):
            batchIds = newIds[i:i + batchSize]
            if len(batchIds) == batchSize or self.isStreamClosed:
                stepId = self._insertFunctionStep(self.samplingStep, batchIds, needsGPU=False,
                                              prerequisites=[])
                self.insertedIds.update(batchIds)
                deps.append(stepId)

        return deps

    def _restoreRuntimeStateFromFinishedSteps(self, doneIds):
        doneIds = set(doneIds)
        self.sampleIds.update(doneIds)
        self.insertedIds.update(doneIds)
        self.processedIds.update(doneIds)

        for step in self._steps:
            isFinished = getattr(step, "isFinished", None)
            if not callable(isFinished) or not isFinished():
                continue

            funcName = getattr(step, "funcName", None)
            if callable(getattr(funcName, "get", None)):
                funcName = funcName.get()

            if funcName != "samplingStep":
                continue

            processedIds = set()
            sampledIds = set()

            argsStr = getattr(step, "argsStr", None)
            if callable(getattr(argsStr, "get", None)):
                argsStr = argsStr.get()

            if argsStr:
                try:
                    stepArgs = json.loads(argsStr)
                    if stepArgs:
                        processedIds.update(stepArgs[0])
                except (TypeError, ValueError):
                    pass

            resultFiles = getattr(step, "_resultFiles", None)
            if callable(getattr(resultFiles, "hasValue", None)):
                if not resultFiles.hasValue():
                    resultFiles = None

            if callable(getattr(resultFiles, "get", None)):
                resultFiles = resultFiles.get()

            if resultFiles:
                try:
                    resultFiles = json.loads(resultFiles)
                except (TypeError, ValueError):
                    resultFiles = []

                for stateFile in resultFiles:
                    if not stateFile or not os.path.exists(stateFile):
                        continue

                    try:
                        with open(stateFile, "r", encoding="utf-8") as handle:
                            state = json.load(handle)
                    except (OSError, TypeError, ValueError):
                        continue

                    processedIds.update(state.get("processedIds", []))
                    sampledIds.update(state.get("sampledIds", []))
                    break

            # Old runs did not persist the pending random selection.
            # Samples already published are still recoverable from the output.
            sampledIds.update(processedIds.intersection(doneIds))

            self.insertedIds.update(processedIds)
            self.processedIds.update(processedIds)
            self.sampleIds.update(sampledIds)

    def _getSamplingStateFile(self, newIds):
        batchKey = hashlib.sha1(
            json.dumps(
                list(newIds),
                separators=(",", ":"),
            ).encode("utf-8")
        ).hexdigest()

        return self._getExtraPath(
            "data_sampler_%s.json" % batchKey
        )

    def samplingStep(self, newIds):
        proportion = self.samplingProportion.get()
        sampledIds = sample_proportion(newIds, proportion)
        stateFile = self._getSamplingStateFile(newIds)
        tmpStateFile = stateFile + ".tmp"

        os.makedirs(
            os.path.dirname(stateFile),
            exist_ok=True,
        )

        state = {
            "processedIds": sorted(newIds),
            "sampledIds": sorted(sampledIds),
        }

        try:
            with open(tmpStateFile, "w", encoding="utf-8") as handle:
                json.dump(state, handle)
            os.replace(tmpStateFile, stateFile)
        finally:
            if os.path.exists(tmpStateFile):
                os.remove(tmpStateFile)

        self.processedIds.update(newIds)
        self.sampleIds.update(sampledIds)

        self.info('From %d new images, %d were random sampled with a proportion of %.2f'
                  % (len(newIds), len(sampledIds), proportion))

        return stateFile


    # ------------------------- UTILS functions --------------------------------
    def _getAllDoneIds(self):
        doneIds = []
        sizeOutput = 0

        if hasattr(self, OUTPUT):
            sizeOutput = self.outputSet.getSize()
            doneIds.extend(list(self.outputSet.getIdSet()))

        return doneIds, sizeOutput

    def _summary(self):
        pass


def sample_proportion(ids_list, proportion=0.25):
    """
    Randomly sample a proportion of elements from a list.

    Parameters:
    - ids_list (list): The list of elements to sample from.
    - proportion (float): The proportion of the list to sample (e.g., 0.25 for 25%).

    Returns:
    - list: A new list containing the sampled elements.
    """
    # Calculate the number of elements to sample
    sample_size = int(len(ids_list) * proportion)

    # Randomly sample elements without replacement
    sampled_list = random.sample(ids_list, sample_size)

    return sampled_list
