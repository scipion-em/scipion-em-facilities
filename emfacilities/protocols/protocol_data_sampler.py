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
import random

from pyworkflow import VERSION_3_0
from pwem.objects import SetOfImages, Set
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils

from pwem.protocols import EMProtocol
from pyworkflow import UPDATED, NEW
from pyworkflow.protocol.constants import STATUS_NEW



OUTPUT = "outputSet"

class ProtDataSampler(EMProtocol):
    """
    Creates randomized subsets from streaming or static image datasets in
    order to reduce data volume while preserving representative sampling
    across large cryo-EM workflows.

    AI Generated:

    Data Sampler (ProtDataSampler) - User Manual
        Overview

        The Data Sampler protocol is designed to generate smaller,
        randomly selected subsets from large image collections such as
        particles, micrographs, or other image-based cryo-EM datasets.
        Its main purpose is to reduce dataset size in a controlled and
        statistically representative way while preserving enough data for
        exploratory analysis, rapid testing, benchmarking, or preliminary
        processing workflows.

        In practical cryo-EM environments, datasets can become extremely
        large, especially during streaming acquisition or automated
        processing sessions. Working with the complete dataset may be
        computationally expensive or unnecessary during early evaluation
        stages. This protocol allows users to extract a fraction of the
        available images while maintaining random selection across the
        incoming data population.

        General Workflow

        The protocol receives an input image collection and continuously
        monitors it during execution. Images are grouped into batches,
        and once enough new data becomes available, a random subset is
        selected according to the user-defined sampling proportion. The
        selected images are then transferred to the output dataset for
        downstream processing.

        This strategy is particularly useful in streaming workflows where
        data acquisition is still ongoing. Instead of waiting for the
        entire experiment to finish, users can immediately begin testing
        classification, reconstruction, or quality-control procedures on
        representative subsets.

        Input Data Considerations

        The protocol accepts any image-based dataset derived from the
        Scipion image framework, including particles, micrographs, and
        related image collections. Since the sampling process is random,
        the biological interpretation of the resulting subset depends on
        the diversity and quality of the original dataset.

        For highly heterogeneous samples, random sampling generally
        preserves the overall population distribution when enough images
        are selected. However, very small sampling proportions may fail
        to capture rare conformational states, uncommon particle views,
        or low-abundance structural populations. Users interested in
        detecting subtle heterogeneity should therefore choose sampling
        proportions carefully.

        Batch Size and Streaming Behavior

        The batch size determines how many new images must accumulate
        before sampling is performed. Smaller batch sizes allow more
        responsive streaming behavior and faster early feedback during
        data collection, but they may increase processing overhead and
        produce noisier statistical representation.

        Larger batch sizes improve statistical stability because the
        random selection occurs over a broader image population. This is
        often preferable for large-scale production workflows or when the
        sampled subset will be used for biologically meaningful
        interpretation.

        During live acquisition, the protocol continuously evaluates the
        incoming dataset and processes only newly available images. This
        allows long-running experiments to generate progressively updated
        sampled outputs without reprocessing previously handled data.

        Sampling Proportion

        The sampling proportion controls how much of each batch is kept.
        A value close to one preserves most of the dataset, while smaller
        values aggressively reduce dataset size. The optimal value depends
        on the intended downstream application.

        For rapid testing of processing parameters, very small subsets
        are often sufficient and can dramatically reduce computational
        cost. For structural interpretation or classification tasks,
        larger sampling proportions are generally recommended to preserve
        biological diversity and angular coverage.

        Random selection is especially useful for creating unbiased test
        datasets, validating workflows, benchmarking algorithms, or
        performing quick quality assessments during microscope sessions.

        Outputs and Interpretation

        The protocol produces a new image dataset containing only the
        randomly selected subset. The output preserves the metadata and
        structural organization required for downstream cryo-EM
        processing pipelines.

        Because the selection is random, repeated executions may produce
        different subsets even when applied to the same dataset. This is
        biologically acceptable in most exploratory workflows, although
        users performing strict reproducibility studies may wish to
        control randomness externally.

        In streaming conditions, the output dataset grows progressively
        over time until the input stream is closed and all eligible
        images have been evaluated.

        Practical Recommendations

        For rapid workflow validation, users commonly begin with small
        sampling proportions and moderate batch sizes. This provides fast
        turnaround while still preserving enough diversity for testing
        alignment, classification, or reconstruction parameters.

        For heterogeneous samples or difficult datasets, larger sampling
        proportions are advisable to avoid unintentionally excluding rare
        structural states. When downstream analysis depends strongly on
        particle diversity, users should visually inspect the sampled
        dataset before drawing biological conclusions.

        In facility or automated processing environments, this protocol
        can substantially reduce computational cost by limiting the
        number of images entering expensive downstream steps during early
        exploratory analysis.

        Final Perspective

        Randomized dataset reduction is an important strategy in modern
        cryo-EM processing because it allows efficient exploration of
        large datasets without requiring full-scale computation at every
        stage. By producing representative subsets during streaming or
        offline processing, the Data Sampler protocol enables faster
        experimentation, rapid quality assessment, and more efficient
        allocation of computational resources while preserving the
        biological relevance of the sampled data.
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
        # Check if there are new images to process from the input set
        self.lastCheck = getattr(self, 'lastCheck', datetime.now())
        mTime = datetime.fromtimestamp(os.path.getmtime(self.inputFn))
        self.debug('Last check: %s, modification: %s'
                    % (pwutils.prettyTime(self.lastCheck),
                        pwutils.prettyTime(mTime)))
        # If the input.sqlite have not changed since our last check,
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
            doneIdsSet = set(doneIds)
            newIdsSet = set(newIds)
            skipIds = list(newIdsSet & doneIdsSet)
            newIds = list(newIdsSet - doneIdsSet)
            self.info("Skipping Images with ID: %s, seems to be done" % skipIds)
            self.insertedIds = set(doneIds)  # During the first round of "Continue" action it has to be filled

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
        doneIdSet = set(doneListIds)
        newDone = list(self.sampleIds - doneIdSet)
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
                self.insertedIds.update(batchIds)
                deps.append(stepId)

        return deps

    def samplingStep(self, newIds):
        proportion = self.samplingProportion.get()
        sampledIds = sample_proportion(newIds, proportion)
        self.processedIds.update(newIds)
        self.sampleIds.update(sampledIds)

        self.info('From %d new images, %d were random sampled with a proportion of %.2f'
                  %(len(newIds), len(sampledIds), proportion))

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