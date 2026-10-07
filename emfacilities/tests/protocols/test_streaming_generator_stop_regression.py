# ***************************************************************************
# * Authors:    Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
# ***************************************************************************/
"""The streaming generator must abandon its loop once the run is over.

A failed step makes pyworkflow mark the protocol as FAILED and the step
executor break out of its own loop - and then join every thread it
started, the steps generator among them. A generator that keeps polling
is never joined, so the run hangs for good with nothing left to do.
"""
import threading
import unittest

import pyworkflow.protocol.constants as cons
from pyworkflow.object import Set

from emfacilities.protocols.protocol_streaming_base import (
    ProtFacilitiesStreamingBase,
)
from emfacilities.protocols.protocol_good_classes_extractor import (
    ProtGoodClassesExtractor,
)


# A real hang cannot be asserted on, so every fake caps its own polling:
# the cap being reached is the failure signal.
MAX_POLLS = 50


class _Status:
    """Mimics the String attribute pyworkflow keeps the status in."""

    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value


class _LoopHarness:
    """Minimal stand-in for a protocol running the shared streaming loop."""

    # The guard under test, taken from the real base class rather than
    # reimplemented here.
    _streamingMustStop = ProtFacilitiesStreamingBase._streamingMustStop

    def __init__(self, status=cons.STATUS_RUNNING, finishAfter=None):
        self.status = _Status(status)
        self.finished = False
        self.polls = 0
        self.closedOutput = 0
        self.sleeps = 0
        self._finishAfter = finishAfter

    def _onStreamingIteration(self):
        self.polls += 1
        if self.polls > MAX_POLLS:
            raise AssertionError(
                "The generator polled %d times without leaving its loop: "
                "the executor's join() would never return." % self.polls
            )

    def _checkNewInput(self):
        pass

    def _checkNewOutput(self):
        if self._finishAfter is not None and self.polls >= self._finishAfter:
            self.finished = True

    def _streamingSleepOnWait(self):
        self.sleeps += 1

    def _closeOutputSet(self):
        self.closedOutput += 1


class TestSharedLoopStopsWhenTheRunIsOver(unittest.TestCase):
    """_runStreamingLoop backs the data counter and the data sampler."""

    def _runLoop(self, harness):
        ProtFacilitiesStreamingBase._runStreamingLoop(harness)

    def testFailedRunLeavesTheLoopImmediately(self):
        harness = _LoopHarness(status=cons.STATUS_FAILED)

        self._runLoop(harness)

        self.assertEqual(
            harness.polls,
            0,
            "A FAILED protocol must not poll even once: the executor is "
            "already joining this very thread.",
        )

    def testAbortedRunLeavesTheLoopImmediately(self):
        harness = _LoopHarness(status=cons.STATUS_ABORTED)

        self._runLoop(harness)

        self.assertEqual(harness.polls, 0)

    def testFailureMidStreamStopsTheLoop(self):
        harness = _LoopHarness(status=cons.STATUS_RUNNING)
        realCheck = harness._checkNewOutput

        def failOnThirdPoll():
            realCheck()
            if harness.polls == 3:
                harness.status = _Status(cons.STATUS_FAILED)

        harness._checkNewOutput = failOnThirdPoll

        self._runLoop(harness)

        self.assertEqual(
            harness.polls,
            3,
            "The loop must stop on the poll right after the failure, not "
            "keep spinning until the stream closes on its own.",
        )

    def testHealthyRunStillPollsUntilItFinishes(self):
        harness = _LoopHarness(status=cons.STATUS_RUNNING, finishAfter=4)

        self._runLoop(harness)

        self.assertEqual(
            harness.polls,
            4,
            "A running protocol must keep polling: the stop guard must "
            "not short-circuit a healthy stream.",
        )
        self.assertEqual(harness.closedOutput, 1)

    def testAPlainStringStatusIsUnderstood(self):
        """status is not always wrapped: accept the bare value too."""
        harness = _LoopHarness()
        harness.status = cons.STATUS_FAILED

        self._runLoop(harness)

        self.assertEqual(harness.polls, 0)

    def testAMissingStatusDoesNotBreakTheLoop(self):
        harness = _LoopHarness(finishAfter=1)
        del harness.status

        self._runLoop(harness)

        self.assertEqual(harness.polls, 1)


class _ExtractorHarness:
    """Minimal stand-in for the good classes extractor generator."""

    LIST_CLASSES = ProtGoodClassesExtractor.LIST_CLASSES
    LIST_IDS = ProtGoodClassesExtractor.LIST_IDS
    _streamingMustStop = ProtFacilitiesStreamingBase._streamingMustStop

    def __init__(self, status=cons.STATUS_RUNNING):
        self.status = _Status(status)
        self.finished = False
        self.isStreamClosed = Set.STREAM_OPEN
        self.newDeps = []
        self._lock = threading.Lock()
        self.polls = 0
        self.insertedFuncs = []
        self.initialised = 0

    # -- collaborators the generator leans on -----------------------------
    def initialStep(self):
        self.initialised += 1

    def _newParticlesToProcess(self):
        self.polls += 1
        if self.polls > MAX_POLLS:
            raise AssertionError(
                "The generator polled %d times without leaving its loop: "
                "the executor's join() would never return." % self.polls
            )
        return True

    def _loadInputClassesSet(self):
        return _ClassSetStub(self.isStreamClosed)

    def _insertFunctionStep(self, func, *args, **kwargs):
        self.insertedFuncs.append(getattr(func, '__name__', str(func)))
        return len(self.insertedFuncs)

    def _streamingSleepOnWait(self):
        pass

    def info(self, *args):
        pass

    def selectGoodClasses(self):
        pass

    def extractElements(self, inputClasses):
        pass

    def closeOutputStep(self):
        pass


class _ClassSetStub:
    def __init__(self, streamState):
        self._streamState = streamState

    def getStreamState(self):
        return self._streamState

    def close(self):
        pass


class TestExtractorGeneratorStopsWhenTheRunIsOver(unittest.TestCase):
    """The extractor keeps a polling loop of its own, so it needs the
    very same guard the shared loop got."""

    def _runGenerator(self, harness):
        ProtGoodClassesExtractor.stepsGeneratorStep(harness)

    def testFailedRunLeavesTheLoopImmediately(self):
        harness = _ExtractorHarness(status=cons.STATUS_FAILED)

        self._runGenerator(harness)

        self.assertEqual(
            harness.polls,
            0,
            "A FAILED protocol must not poll even once: the executor is "
            "already joining this very thread.",
        )

    def testAbortedRunLeavesTheLoopImmediately(self):
        harness = _ExtractorHarness(status=cons.STATUS_ABORTED)

        self._runGenerator(harness)

        self.assertEqual(harness.polls, 0)

    def testFailureMidStreamStopsTheLoop(self):
        harness = _ExtractorHarness()
        realCheck = harness._newParticlesToProcess

        def failOnThirdPoll():
            hasNew = realCheck()
            if harness.polls == 3:
                harness.status = _Status(cons.STATUS_FAILED)
            return hasNew

        harness._newParticlesToProcess = failOnThirdPoll

        self._runGenerator(harness)

        self.assertEqual(
            harness.polls,
            3,
            "The loop must stop on the poll right after the failure, not "
            "keep spinning until the input classes stream closes.",
        )

    def testHealthyRunStillReachesTheStreamClose(self):
        harness = _ExtractorHarness()

        def closeOnThirdPoll():
            harness.polls += 1
            if harness.polls >= 3:
                harness.isStreamClosed = Set.STREAM_CLOSED
            return True

        harness._newParticlesToProcess = closeOnThirdPoll

        self._runGenerator(harness)

        self.assertTrue(
            harness.finished,
            "A running protocol must still poll until the producer "
            "closes: the stop guard must not short-circuit it.",
        )
        self.assertIn('closeOutputStep', harness.insertedFuncs)

    def testGoodClassesAreSelectedOnlyOncePerRun(self):
        """Whether the selection ran is answered by another thread.

        Polling that flag means the generator can queue the very same
        selection several times before the first one gets to run, and
        every extraction step would then wait on a redundant step.
        """
        harness = _ExtractorHarness()

        def closeOnFourthPoll():
            harness.polls += 1
            if harness.polls >= 4:
                harness.isStreamClosed = Set.STREAM_CLOSED
            return True

        harness._newParticlesToProcess = closeOnFourthPoll

        self._runGenerator(harness)

        selections = harness.insertedFuncs.count('selectGoodClasses')
        self.assertEqual(
            selections,
            1,
            "The selection step was inserted %d times: the generator must "
            "not wait for another thread to flip the flag." % selections,
        )


if __name__ == '__main__':
    unittest.main()
