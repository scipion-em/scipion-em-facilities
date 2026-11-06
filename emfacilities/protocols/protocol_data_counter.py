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

from .protocol_streaming_base import ProtFacilitiesStreamingBase
from pwem.objects import SetOfImages, Set

from pyworkflow import VERSION_3_0
import pyworkflow.protocol.params as params
from pyworkflow import UPDATED, NEW


OUTPUT = "outputSet"


class ProtDataCounter(ProtFacilitiesStreamingBase):
    """
    Collects and transfers a controlled number of images from streaming
    or static cryo-EM datasets, allowing workflows to stop automatically
    once a desired dataset size or acquisition time limit has been reached.

    AI Generated:

    Data Counter (ProtDataCounter) - User Manual
        Overview

        The Data Counter protocol is designed to monitor incoming image
        datasets and generate an output collection once a predefined
        number of images has been accumulated. Its primary purpose is to
        provide controlled termination points in streaming or automated
        cryo-EM workflows, enabling users to limit dataset size or stop
        acquisition-driven processing after a specified amount of time.

        In practical cryo-EM environments, data acquisition can continue
        for many hours or even days. During this process, users may only
        require a subset of the full dataset for rapid testing,
        preliminary reconstruction, benchmarking, or quality assessment.
        This protocol provides a simple mechanism to automatically close
        processing streams once sufficient data has been collected.

        General Workflow

        The protocol continuously monitors an incoming image collection
        such as particles, micrographs, or related cryo-EM image sets.
        As new images become available, they are progressively registered
        and transferred to the output dataset until the desired stopping
        condition is reached.

        Two operating modes are available. In the standard mode, the
        protocol stops only after the requested number of images has been
        accumulated. In timer-assisted mode, the protocol can also stop
        automatically when a predefined acquisition time limit has
        expired, even if the target number of images has not yet been
        reached.

        This dual behavior is particularly useful in automated facilities
        and live microscopy sessions where users may need to balance data
        quantity, computational cost, and available acquisition time.

        Input Data Considerations

        The protocol accepts any image-based dataset compatible with the
        Scipion image framework, including particles, micrographs, and
        related image collections. Since the protocol transfers incoming
        images progressively, the biological representativeness of the
        output depends on the order and diversity of the acquired data.

        For homogeneous samples, early subsets are often sufficient for
        rapid validation of acquisition quality and processing stability.
        For heterogeneous or conformationally flexible systems, users
        should consider that stopping acquisition too early may exclude
        rare structural populations or poorly represented particle views.

        Output Size Control

        The output size parameter defines the maximum number of images
        that will be transferred to the output dataset. Once this limit
        is reached, the protocol finalizes the stream and closes the
        output collection.

        This functionality is especially valuable when preparing reduced
        datasets for rapid experimentation or when computational
        resources are limited. Users can quickly generate manageable
        subsets suitable for testing alignment, classification, or
        reconstruction strategies without processing the entire dataset.

        In high-throughput environments, limiting dataset size also helps
        reduce storage consumption and accelerates downstream processing.

        Timer-Based Streaming Control

        An optional timer mechanism allows the protocol to terminate
        acquisition-driven workflows after a user-defined period of time.
        This is useful when microscope access is limited, when overnight
        acquisition sessions must stop automatically, or when users wish
        to evaluate early data quality before committing to longer runs.

        The timer accepts flexible human-readable durations including
        seconds, minutes, hours, or days. Once the timer expires, the
        protocol closes the output stream even if the requested image
        count has not been reached.

        From a practical perspective, this mode is valuable for facility
        pipelines where rapid preliminary feedback is more important than
        complete dataset collection.

        Streaming Behavior

        During execution, the protocol continuously checks for new input
        images and progressively updates the output dataset. This allows
        downstream workflows to begin immediately while acquisition is
        still ongoing.

        The protocol is designed to function safely in live streaming
        environments, ensuring that all available data is properly
        registered before the stream is finalized. This behavior is
        particularly important in automated cryo-EM pipelines where data
        may continue arriving asynchronously during acquisition.

        Outputs and Interpretation

        The protocol produces a new image dataset containing the selected
        subset of images accumulated before the stopping condition was
        reached. The output preserves the metadata and organizational
        structure required for downstream cryo-EM analysis.

        When the timer mode is enabled, the final dataset size may differ
        from the requested target if the acquisition period ends before
        enough images are collected. Conversely, when acquisition is
        sufficiently fast, the output size condition may terminate the
        workflow before the timer expires.

        Practical Recommendations

        For rapid quality-control workflows, relatively small output
        sizes are often sufficient to evaluate particle quality,
        microscope alignment, ice thickness, or reconstruction stability.
        This enables fast iteration during microscope sessions without
        waiting for complete acquisition.

        In heterogeneous biological systems, users should avoid selecting
        output sizes that are too small, as this may reduce angular
        coverage or exclude low-abundance conformational states. Larger
        subsets generally provide more biologically reliable
        representation.

        Timer mode is particularly useful for unattended overnight runs,
        automated facility pipelines, and preliminary exploratory
        experiments where acquisition time must remain bounded.

        Final Perspective

        Controlled streaming termination is an important component of
        modern cryo-EM automation because it allows users to balance data
        quantity, computational cost, and experimental timing. By
        automatically stopping workflows according to dataset size or
        acquisition duration, the Data Counter protocol supports more
        efficient exploratory analysis, faster quality assessment, and
        better resource management across large-scale cryo-EM pipelines.
    """

    _label = 'data counter'
    _devStatus = NEW
    _lastUpdateVersion = VERSION_3_0
    _possibleOutputs = {OUTPUT: SetOfImages}

    def __init__(self, **args):
        ProtFacilitiesStreamingBase.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputImages', params.PointerParam, pointerClass='SetOfImages',
                      label="Input images", important=True)
        form.addParam('outputSize', params.IntParam, default=10000,
                      label='Output size',
                      help='How many images need to be on input to '
                           'create output set.')

        form.addSection(label='Stream of data timer')
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

        self._defineStreamingParams(form)
        form.addParallelSection(threads=3, mpi=1)

# --------------------------- INSERT steps functions -------------------------
    def stepsGeneratorStep(self):
        self.initializeParams()
        self._runStreamingLoop()

    def _onStreamingIteration(self):
        if self.boolTimer.get() and not self.timerOut:
            self.timerStep()

    def initializeParams(self):
        self.finished = False
        # Important to have both:
        self.insertedIds = set() # Contains images that have been inserted in a Step (checkNewInput).
        self.processedIds = set() # Ids to be output
        # Discovery watermark only. It is intentionally rebuilt from zero on
        # Continue; persisted outputs remain the source of truth for completion.
        self._lastInputId = 0
        self.isStreamClosed = self.inputImages.get().isStreamClosed()
        # Contains images that have been processed in a Step (checkNewOutput).
        self._initInputTypeState()
        self.limitReach = False
        self.timerOut = False
        self.timeoutSecs = self.getTimeOutInSeconds(self.timeout.get())
        self.lastTimeCheckTimer = datetime.now() # Timer
        self.lastRound = False

    def _checkNewInput(self):
        # Check if there are new images to process from the input set
        if self.finished:
            return

        # Always inspect the logical input set. Backing-file mtimes are not
        # a reliable change detector for streamed logical Set contents.
        if self.lastRound:
            self.info("Last round sleeping for 10 seconds to allow all the input to be loaded")
            time.sleep(10) # Needs to make sure that eventhough the stream is closed all the data in the inputset is loaded

        newIds, producerClosed = self._discoverNewInputIds(self.insertedIds)
        self.lastCheck = datetime.now()
        # Keep the historical terminal wait while the input Set catches
        # up; _discoverNewInputIds already refuses to declare the consumer
        # stream closed until every item advertised by getSize() is visible.
        self.lastRound = producerClosed

        if self.isContinued() and not self.insertedIds:  # For "Continue" action and the first round
            doneIds, _ = self._getAllDoneIds()
            skipIds = list(set(newIds).intersection(set(doneIds)))
            newIds = list(set(newIds).difference(set(doneIds)))
            self.info("Skipping Images with ID: %s, seems to be done" % skipIds)
            self.insertedIds = set(doneIds) # During the first round of "Continue" action it has to be filled
            # Same reconciliation seeds what is already published, so the
            # output poll never has to ask the output itself.
            self._markOutputIdsPersisted(OUTPUT, doneIds)

        if newIds and not self.limitReach and not self.timerOut:
            self._insertNewImageSteps(newIds)

    def _checkNewOutput(self):
        if self.finished:
            return

        # During normal polling processedIds is only the in-memory queue of
        # items whose processing step finished but whose output has not yet
        # been committed. Avoid rebuilding the complete persisted ID set here:
        # Continue performs that reconciliation separately in _checkNewInput.
        currentOutputSize = (
            self.outputSet.getSize() if hasattr(self, OUTPUT) else 0
        )
        newDone = set(self.processedIds)
        allDone = currentOutputSize + len(newDone)
        limitOutputSize = self.outputSize.get()

        inputSet = self._loadInputSet(None)
        try:
            maxSize = inputSet.getSize()
            self.limitReach = allDone >= limitOutputSize

            # We have finished when there is not more input images
            # (stream closed) or when the limit of output size is met
            self.finished = (self.isStreamClosed and allDone == maxSize) or (self.limitReach or self.timerOut)

            if not self.finished and not newDone:
                # If we are not finished and no new output have been produced
                # it does not make sense to proceed and updated the outputs
                # so we exit from the function here
                return

            outputSet = self._loadOutputSet(self._inputClass, self._baseName,
                                            outputName=OUTPUT)

            if currentOutputSize < limitOutputSize:
                # A commit can write and still be reported as failed,
                # and processedIds is only emptied once it is known to
                # have succeeded - so an id can come round again after
                # it is already in the output. Appending it twice would
                # duplicate it.
                alreadyPublished = self._getKnownPersistedOutputIds(OUTPUT)

                persistedNow = set()
                for imageId in newDone:
                    if imageId in alreadyPublished:
                        persistedNow.add(imageId)
                        continue

                    # Set.getItem raises rather than returning None for a
                    # row it cannot find. An id discovered earlier via
                    # _discoverIdsAfter is not guaranteed to still be
                    # selectable on this freshly reloaded input Set: a
                    # row can lag behind the metadata that advertised it.
                    # Check membership first and leave it pending for the
                    # next round instead of crashing the whole protocol.
                    if imageId not in inputSet:
                        self.error(
                            "Image with id %d is not yet visible in the "
                            "input Set; leaving it pending for the next "
                            "round." % imageId
                        )
                        continue

                    image = inputSet.getItem("id", imageId).clone()
                    outputSet.append(image)
                    persistedNow.add(imageId)
                    currentOutputSize += 1
                    if currentOutputSize == limitOutputSize:
                        break # We have reach the limit for the outputSize

                # Recompute completion from what was ACTUALLY persisted this
                # round, not the optimistic pre-loop count: an item that is
                # not yet visible must not let the protocol close the
                # output before it is actually persisted.
                self.limitReach = currentOutputSize >= limitOutputSize
                self.finished = (
                    (self.isStreamClosed and currentOutputSize == maxSize)
                    or self.limitReach
                    or self.timerOut
                )

                streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

                self._updateOutputSet(OUTPUT, outputSet, streamMode)
                # Only forget pending IDs after the output update succeeds.
                # If persistence raises, they remain queued for retry.
                self._markOutputIdsPersisted(OUTPUT, persistedNow)
                self.processedIds.difference_update(persistedNow)
        finally:
            inputSet.close()

        self._store()

    def _insertNewImageSteps(self, newIds):
        """ Insert steps to register new images (from streaming)
        Params:
            newIds: input images ids to be processed
        """
        deps = []
        stepId = self._insertFunctionStep(self.registerStep, newIds, needsGPU=False,
                                          prerequisites=[])
        deps.append(stepId)
        self.insertedIds.update(newIds)

        return deps

    def registerStep(self, newIds):
        self.info('Registering the %d new images' % len(newIds))
        self.processedIds.update(newIds)


    def timerStep(self):
        now = datetime.now()

        if self.initTime.hasValue():
            startTime = self.initTime.datetime()
            timeoutSecs = self.getTimeOutInSeconds(self.timeout.get())
            endTime = startTime + timedelta(seconds=timeoutSecs)
        else:
            # Fallback for isolated/unit usage where the protocol has not
            # gone through Protocol.setRunning().
            endTime = self.lastTimeCheckTimer + timedelta(seconds=self.timeoutSecs)

        remainingTime = (endTime - now).total_seconds()

        if remainingTime <= 0:
            self.timeoutSecs = 0
            self.timerOut = True
            self.info("  timer is consumed terminating protocol.")
            self.summaryVar.set("Timer is consumed terminating protocol.")
        else:
            self.timeoutSecs = int(remainingTime)
            self.info(f"  remaining time: {self.timeoutSecs} seconds.")
            self.summaryVar.set(
                "Time activated remaining time: %d seconds" % self.timeoutSecs
            )

        self.lastTimeCheckTimer = now


    # ------------------------- UTILS functions --------------------------------
    def getTimeOutInSeconds(self, timeOut):
        timeOutFormatRegexList = {r'\d+s': 1, r'\d+m': 60, r'\d+h': 3600,
                                  r'\d+d': 86400}
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
        summary = []
        if not hasattr(self, OUTPUT):
            summary.append("Output set not ready yet.")
        else:
            outputExpectedSize = self.outputSize.get()
            outputSize = self.outputSet.getSize()
            per = (outputSize/outputExpectedSize) * 100

            summary.append("%.1f %% of the expected output length: %d / %d "
                           % (per, outputSize, outputExpectedSize))

        if self.boolTimer.get():
            summary.append(self.summaryVar.get())

        return summary
