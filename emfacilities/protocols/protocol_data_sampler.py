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
from datetime import datetime
import time
import copy
import numpy as np
import random
from collections import defaultdict

from pyworkflow import VERSION_3_0
from pwem.objects import SetOfImages, Set
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils

from pwem.protocols import EMProtocol
from pyworkflow import UPDATED, NEW
from pyworkflow.protocol.constants import STATUS_NEW, STEPS_PARALLEL



OUTPUT = "outputSet"

class ProtDataSampler(EMProtocol):
    """
    Protocol to make a subset of images from the original one.
    Waits until certain batch of images is prepared and then it samples a percentage of it and send them to output
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
                           'the random sampling.')
        form.addSection(label='Stream of data timer')
        form.addParam('delay', params.IntParam, default=10, label="Delay (sec)",
                      validators=[params.GT(2, "must be larger than 3sec.")],
                      help="Delay in seconds before checking new output")


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
        self.processedIds = [] # Ids to be register to output
        self.sampleIds = [] # Ids to be output
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
        # Check if there are new ctf to process from the input set
        self.lastCheck = getattr(self, 'lastCheck', datetime.now())
        mTime = datetime.fromtimestamp(os.path.getmtime(self.inputFn))
        self.debug('Last check: %s, modification: %s'
                    % (pwutils.prettyTime(self.lastCheck),
                        pwutils.prettyTime(mTime)))
        # If the input movies.sqlite have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime and self.insertedIds:  # If this is empty it is dut to a static "continue" action or it is the first round
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

        # Now handle the steps depending on the streaming batch size
        batchSize = self.batchSize.get()
        if len(newIds) < batchSize and not self.isStreamClosed:
            return  # No register any step if the batch size is not reach unless is the lass iter

        if newIds:
            fDeps = self._insertNewImageSteps(newIds, batchSize)
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()

    def _checkNewOutput(self):
        doneListIds, currentOutputSize = self._getAllDoneIds()
        #processedIds = copy.deepcopy(self.processedIds)
        sampleIds = copy.deepcopy(self.sampleIds)
        newDone = [imageId for imageId in sampleIds if imageId not in doneListIds]
        allDone = len(doneListIds) + len(newDone)
        maxSize = int(self._loadInputSet(self.inputFn).getSize() * self.samplingProportion.get())

        # We have finished when there is not more input images
        # (stream closed) or when the limit of output size is met
        self.finished = self.isStreamClosed and allDone == maxSize
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if not self.finished and not newDone:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        inputSet = self._loadInputSet(self.inputFn)
        outputSet = self._loadOutputSet(self._inputClass, self._baseName)

        for imageId in newDone:
            image = inputSet.getItem("id", imageId).clone()
            outputSet.append(image)

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
                self.insertedIds.extend(batchIds)
                deps.append(stepId)

        return deps

    def samplingStep(self, newIds):
        proportion = self.samplingProportion.get()
        sampledIds = sample_proportion(newIds, proportion)
        self.processedIds.extend(newIds)
        self.sampleIds.extend(sampledIds)

        self.info('From %d new images, %d were random sampled with a proportion of %2f.'
                  %(len(newIds), len(sampledIds), proportion))

        if not self.isStreamClosed:
            self.delayRegister()

    def delayRegister(self):
        delay = self.delay.get()
        self.info('Sleeping for delay time: %d' %delay)
        time.sleep(delay)

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