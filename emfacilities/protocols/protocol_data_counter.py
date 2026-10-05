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
    Protocol to make a subset of images from the original one. Waits until certain number of images is prepared and then send them to output.
    It can works in 2 ways:
        - Simple mode: once the number of items is reached, a setOfImages is returned and
            the protocol finishes (ending the streaming from this point).
        - If timer activated: either once the number of items is reached, or the timer is consumed
            a setOfImages is returned and the protocol finishes (ending the streaming from this point).
    The protocol will accept Micrographs, Particles and any kind of object that inherits from Image base class.
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

        while not self.finished:
            if self.boolTimer.get() and not self.timerOut:
                self.timerStep()

            self._checkNewInput()
            self._checkNewOutput()

            if not self.finished:
                self._streamingSleepOnWait()

        self._closeOutputSet()

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
        self._inputClass = self.inputImages.get().getClass()
        self._inputType = self.inputImages.get().getClassName().split('SetOf')[1]
        self._baseName = '%s.sqlite' % self._inputType.lower()
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

        inputSet = self._loadInputSet(None)
        try:
            newIds, self._lastInputId = self._discoverIdsAfter(
                inputSet,
                self._lastInputId,
            )

            self.lastCheck = datetime.now()
            producerClosed = inputSet.isStreamClosed()

            newIds, terminalConsistent = (
                self._reconcileClosedStreamIds(
                    inputSet,
                    newIds,
                    self.insertedIds,
                    producerClosed,
                )
            )

            # Keep the historical terminal wait while PostgreSQL catches up,
            # but do not declare the consumer stream closed until every item
            # advertised by getSize() is actually visible.
            self.lastRound = producerClosed
            self.isStreamClosed = (
                producerClosed
                and terminalConsistent
            )
        finally:
            inputSet.close()

        if self.isContinued() and not self.insertedIds:  # For "Continue" action and the first round
            doneIds, _ = self._getAllDoneIds()
            skipIds = list(set(newIds).intersection(set(doneIds)))
            newIds = list(set(newIds).difference(set(doneIds)))
            self.info("Skipping Images with ID: %s, seems to be done" % skipIds)
            self.insertedIds = set(doneIds) # During the first round of "Continue" action it has to be filled

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
                persistedNow = set()
                for imageId in newDone:
                    # Set.getItem raises rather than returning None for a
                    # row it cannot find. An id discovered earlier via
                    # _discoverIdsAfter is not guaranteed to still be
                    # selectable on this freshly reloaded input Set (e.g.
                    # replication lag under a PostgreSQL-backed
                    # compatibility bridge) - check membership first and
                    # leave it pending for the next round instead of
                    # crashing the whole protocol.
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
