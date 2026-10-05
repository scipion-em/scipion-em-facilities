# ***************************************************************************
# * Authors:    Daniel Marchán (da.marchan@cnb.csic.es)
# *
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
from datetime import datetime, timedelta
from unittest.mock import patch
from pyworkflow.tests import BaseTest, DataSet
from pwem.protocols.protocol_import import ProtImportMicrographs
from pyworkflow.object import Pointer
import pyworkflow.tests as tests
from emfacilities.protocols.protocol_data_counter import ProtDataCounter, OUTPUT


class TestDataCounter(BaseTest):
    """ Test data counter protocol """

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.protImport = cls.runImportMicrographs(cls.micsFn)


    @classmethod
    def runImportMicrographs(cls, micsFn):
        """
        Import Micrographs
        """
        protImport = cls.newProtocol(ProtImportMicrographs,
                                      filesPath=micsFn,
                                      samplingRate=1.237,
                                      voltage=300)
        cls.launchProtocol(protImport)

        return protImport






    def testDetectsNewInputWhenMtimeDoesNotChange(self):
        class InputSet:
            def getUniqueValues(self, attributes, where=None):
                self.assertQuery = (attributes, where)
                return [4, 5, 6]

            def isStreamClosed(self):
                return False

            def close(self):
                pass

        insertedBatches = []

        prot = self.newProtocol(
            ProtDataCounter,
            outputSize=100,
            boolTimer=False,
        )
        prot.finished = False
        prot.inputFn = "input.sqlite"
        prot.insertedIds = {1, 2, 3}
        prot.processedIds = {1, 2, 3}
        prot._lastInputId = 3
        prot.isStreamClosed = False
        prot.lastRound = False
        prot.limitReach = False
        prot.timerOut = False
        prot.isContinued = lambda: False
        prot._loadInputSet = lambda _: InputSet()
        prot._getFirstJoinStep = lambda: None
        prot.updateSteps = lambda: None

        def insertNewImageSteps(newIds):
            insertedBatches.append(list(newIds))
            return []

        prot._insertNewImageSteps = insertNewImageSteps

        with patch(
                "emfacilities.protocols.protocol_data_counter.os.path.getmtime",
                return_value=0,
        ):
            prot._checkNewInput()

        self.assertEqual(
            insertedBatches,
            [[4, 5, 6]],
            "New logical Set items must be detected even when the sqlite "
            "mtime does not change.",
        )


    def testTimeoutDaysUse24Hours(self):
        prot = self.newProtocol(ProtDataCounter)

        self.assertEqual(prot.getTimeOutInSeconds("1d"), 86400)
        self.assertEqual(
            prot.getTimeOutInSeconds("1d 2h 20m 15s"),
            86400 + 2 * 3600 + 20 * 60 + 15,
        )


    def testTimerUsesPreservedProtocolStartOnContinue(self):
        class SummaryVar:
            def __init__(self):
                self.value = None

            def set(self, value):
                self.value = value

        prot = self.newProtocol(
            ProtDataCounter,
            outputSize=100,
            boolTimer=True,
            timeout="10s",
        )
        prot.finished = False
        prot.timerOut = False
        prot.timeoutSecs = 10
        prot.lastTimeCheckTimer = datetime.now()
        prot.summaryVar = SummaryVar()
        prot.initTime.set(datetime.now() - timedelta(seconds=7))

        prot.timerStep()

        self.assertFalse(prot.timerOut)
        self.assertLessEqual(
            prot.timeoutSecs,
            3,
            "Continue must preserve the elapsed timer budget from the original run.",
        )


    def testTimerExpiresWithoutNewInput(self):
        class SummaryVar:
            def __init__(self):
                self.value = None

            def set(self, value):
                self.value = value

        prot = self.newProtocol(
            ProtDataCounter,
            outputSize=100,
            boolTimer=True,
            timeout="10s",
        )
        prot.finished = False
        prot.timerOut = False
        prot.timeoutSecs = 10
        prot.lastTimeCheckTimer = datetime.now() - timedelta(seconds=11)
        prot.summaryVar = SummaryVar()

        # Simulate an idle generator round. The modern streaming scheduler
        # executes the timer from stepsGeneratorStep(), not from _stepsCheck().
        # This test prepares the runtime state manually. Do not let the
        # generator reinitialize it or require a real input pointer.
        prot.initializeParams = lambda: None
        prot._checkNewInput = lambda: None

        def checkNewOutput():
            if prot.timerOut:
                prot.finished = True

        prot._checkNewOutput = checkNewOutput
        prot._streamingSleepOnWait = lambda: None
        prot._closeOutputSet = lambda: None

        prot.stepsGeneratorStep()

        self.assertTrue(
            prot.timerOut,
            "The timer must expire even when no new input batch arrives.",
        )


    def testDataCounter2000(self):
        prot = self._runDataCounter("Counter images till 1", outputSize=1)
        self.assertSetSize(prot.outputSet, size=1)


    def testDataCounter4000(self):
        prot = self._runDataCounter("Counter images till 2", outputSize=2)
        self.assertSetSize(prot.outputSet, size=2)

    def _runDataCounter(cls, label, outputSize):
        protDataSampler = cls.newProtocol(ProtDataCounter,
                                          outputSize=outputSize,
                                          delay=3)
        protDataSampler.inputImages = Pointer(cls.protImport, extended='outputMicrographs')
        protDataSampler.setObjLabel(label)
        cls.launchProtocol(protDataSampler)

        return protDataSampler


class TestDataCounterLoadOutputSet(tests.unittest.TestCase):
    """Lightweight regression tests that need no real project/dataset."""

    def testLoadOutputSetReusesLogicalOutputWithoutBackingFile(self):
        # Regression test: an output that Scipion already knows about
        # (protocol.outputSet) must be reused even when its backing file
        # was never materialized on disk yet. Falling through to "no
        # backing file -> build a fresh, empty Set" would silently discard
        # whatever was already appended to the real logical output.
        prot = ProtDataCounter()

        class ExistingOutputSet:
            def __init__(self):
                self.loadAllPropertiesCalls = 0
                self.enableAppendCalls = 0
                self.copiedFrom = None

            def loadAllProperties(self):
                self.loadAllPropertiesCalls += 1

            def enableAppend(self):
                if not self.loadAllPropertiesCalls:
                    raise AssertionError("Persisted output must be refreshed before enableAppend().")
                self.enableAppendCalls += 1

            def copyInfo(self, inputs):
                self.copiedFrom = inputs

        existingOutputSet = ExistingOutputSet()
        prot.outputSet = existingOutputSet

        with patch(
                "emfacilities.protocols.protocol_data_counter.os.path.exists",
                return_value=False,
        ):
            outputSet = prot._loadOutputSet(
                object, "images.sqlite", outputName=OUTPUT
            )

        self.assertIs(existingOutputSet, outputSet)
        self.assertEqual(1, existingOutputSet.loadAllPropertiesCalls)
        self.assertEqual(1, existingOutputSet.enableAppendCalls)

class TestDataCounterInputSetLifecycleRegression(tests.unittest.TestCase):

    def testOutputPollingDoesNotRescanAllPersistedIds(self):
        class _Value:
            def get(self):
                return 100

        class _Image:
            def __init__(self, objId):
                self.objId = objId

            def clone(self):
                return _Image(self.objId)

        class _InputSet:
            def __init__(self):
                self.closed = False

            def getSize(self):
                return 5

            def __contains__(self, objId):
                return True

            def getItem(self, field, value):
                return _Image(value)

            def close(self):
                self.closed = True

        class _OutputSet:
            STREAM_OPEN = 1
            STREAM_CLOSED = 2

            def __init__(self):
                self.ids = [1, 2, 3]

            def getSize(self):
                return len(self.ids)

            def getIdSet(self):
                raise AssertionError(
                    "_checkNewOutput must not rescan every persisted output ID "
                    "on each streaming poll."
                )

            def append(self, image):
                self.ids.append(image.objId)

        class _Harness:
            finished = False
            processedIds = {4, 5}
            isStreamClosed = False
            timerOut = False
            outputSize = _Value()
            _inputClass = object
            _baseName = "images.sqlite"

            def __init__(self):
                self.inputSet = _InputSet()
                self.outputSet = _OutputSet()
                self.updated = False

            def _getAllDoneIds(self):
                return ProtDataCounter._getAllDoneIds(self)

            def _loadInputSet(self, _):
                return self.inputSet

            def _loadOutputSet(self, SetClass, baseName, outputName=None):
                return self.outputSet

            def _updateOutputSet(self, outputName, outputSet, streamMode):
                self.updated = True

            def _getFirstJoinStep(self):
                return None

            def _store(self):
                pass

        protocol = _Harness()

        ProtDataCounter._checkNewOutput(protocol)

        self.assertTrue(protocol.updated)
        self.assertEqual([1, 2, 3, 4, 5], protocol.outputSet.ids)
        self.assertEqual(set(), protocol.processedIds)
        self.assertTrue(protocol.inputSet.closed)

    def testCheckNewOutputClosesInputSetWhenThereIsNoNewOutput(self):
        class _Value:
            def get(self):
                return 100

        class _InputSet:
            def __init__(self):
                self.closed = False

            def getSize(self):
                return 3

            def close(self):
                self.closed = True

        class _Harness:
            finished = False
            processedIds = set()
            isStreamClosed = False
            timerOut = False
            outputSize = _Value()

            def __init__(self):
                self.inputFn = "input.sqlite"
                self.inputSet = _InputSet()

            def _getAllDoneIds(self):
                return [], 0

            def _loadInputSet(self, inputFn):
                return self.inputSet

        protocol = _Harness()

        ProtDataCounter._checkNewOutput(protocol)

        self.assertTrue(
            protocol.inputSet.closed,
            "The input Set opened by _checkNewOutput must be closed "
            "even when there is no new output to publish.",
        )


class TestDataCounterBackendIndependence(tests.unittest.TestCase):
    def testInputSetAccessUsesLogicalPointerWithoutBackingFilename(self):
        class _Value:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        class _ForbiddenSetConstructor:
            def __init__(self, *args, **kwargs):
                raise AssertionError(
                    "DataCounter must not reconstruct the input Set from a backing filename."
                )

        class _LogicalInputSet:
            def __init__(self):
                self.loadCalls = 0

            def isStreamClosed(self):
                return False

            def getClass(self):
                return _ForbiddenSetConstructor

            def getClassName(self):
                return "SetOfMicrographs"

            def getFileName(self):
                raise AssertionError(
                    "DataCounter must not require a backing filename for the logical input Set."
                )

            def loadAllProperties(self):
                self.loadCalls += 1

        class _Pointer:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        logicalInput = _LogicalInputSet()
        protocol = ProtDataCounter()
        protocol.inputImages = _Pointer(logicalInput)
        protocol.timeout = _Value("1h")

        protocol.initializeParams()

        loadedInput = protocol._loadInputSet("must-not-be-used.sqlite")

        self.assertIs(logicalInput, loadedInput)
        self.assertEqual(1, logicalInput.loadCalls)


class TestDataCounterStreamingScalability(tests.unittest.TestCase):
    def testDiscoveryQueriesOnlyIdsBeyondWatermark(self):
        class _InputSet:
            def __init__(self):
                self.uniqueCalls = []
                self.closed = False

            def getIdSet(self):
                raise AssertionError(
                    "DataCounter must not fetch the complete input ID set on every poll."
                )

            def getUniqueValues(self, attributes, where=None):
                self.uniqueCalls.append((attributes, where))
                return [4, 5, 6]

            def isStreamClosed(self):
                return False

            def close(self):
                self.closed = True

        class _Harness:
            finished = False
            lastRound = False
            isStreamClosed = False
            limitReach = False
            timerOut = False
            insertedIds = {1, 2, 3}
            processedIds = {1, 2, 3}
            _lastInputId = 3

            def __init__(self):
                self.inputSet = _InputSet()
                self.insertedBatches = []

            def _loadInputSet(self, _):
                return self.inputSet

            def _discoverIdsAfter(self, inputSet, lastId):
                return ProtDataCounter._discoverIdsAfter(self, inputSet, lastId)

            def _reconcileClosedStreamIds(
                    self,
                    inputSet,
                    discoveredIds,
                    knownIds,
                    producerClosed,
            ):
                return ProtDataCounter._reconcileClosedStreamIds(
                    self,
                    inputSet,
                    discoveredIds,
                    knownIds,
                    producerClosed,
                )

            def _getFirstJoinStep(self):
                return None

            def isContinued(self):
                return False

            def _insertNewImageSteps(self, newIds):
                self.insertedBatches.append(list(newIds))
                self.insertedIds.update(newIds)
                return []

            def updateSteps(self):
                pass

        protocol = _Harness()

        ProtDataCounter._checkNewInput(protocol)

        self.assertEqual([("id", "id > 3")], protocol.inputSet.uniqueCalls)
        self.assertEqual([[4, 5, 6]], protocol.insertedBatches)
        self.assertEqual(6, protocol._lastInputId)
        self.assertTrue(protocol.inputSet.closed)

class TestDataCounterStreamingArchitecture(tests.unittest.TestCase):
    def testUsesFacilitiesStreamingBaseAndCoreGeneratorInsertion(self):
        from pyworkflow.protocol import ProtStreamingBase, STEPS_PARALLEL
        from emfacilities.protocols.protocol_streaming_base import (
            ProtFacilitiesStreamingBase,
        )

        self.assertTrue(
            issubclass(ProtDataCounter, ProtFacilitiesStreamingBase)
        )
        self.assertTrue(
            issubclass(ProtFacilitiesStreamingBase, ProtStreamingBase)
        )
        self.assertIs(
            ProtDataCounter._insertAllSteps,
            ProtStreamingBase._insertAllSteps,
        )
        self.assertEqual(
            STEPS_PARALLEL,
            ProtDataCounter.stepsExecutionMode,
        )


    def testGeneratorWaitsForOutputCompletionBeforeClosing(self):
        class _Value:
            def get(self):
                return False

        class _Harness:
            boolTimer = _Value()
            timerOut = False

            def __init__(self):
                self.finished = False
                self.iteration = 0
                self.events = []

            def initializeParams(self):
                self.finished = False
                self.events.append("initialize")

            def _checkNewInput(self):
                self.iteration += 1
                self.events.append("input-%d" % self.iteration)

            def _checkNewOutput(self):
                self.events.append("output-%d" % self.iteration)
                if self.iteration == 2:
                    self.finished = True

            def _streamingSleepOnWait(self):
                self.events.append("sleep-%d" % self.iteration)

            def _closeOutputSet(self):
                self.events.append("close")

        protocol = _Harness()

        ProtDataCounter.stepsGeneratorStep(protocol)

        self.assertEqual(
            [
                "initialize",
                "input-1",
                "output-1",
                "sleep-1",
                "input-2",
                "output-2",
                "close",
            ],
            protocol.events,
        )

class TestDataCounterClosedStreamLateVisibilityRegression(
        tests.unittest.TestCase):
    def testClosedStreamReconcilesIdsBelowWatermark(self):
        class _Value:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        class _InputSet:
            def __init__(self):
                self.uniqueCalls = []
                self.closedCalls = 0

            def getUniqueValues(self, attributes, where=None):
                self.uniqueCalls.append((attributes, where))

                if where is None:
                    return list(range(1, 11))

                if where == "id > 0":
                    return [9, 10]

                if where == "id > 10":
                    return []

                raise AssertionError(
                    "Unexpected discovery query: %r" % (where,)
                )

            def getSize(self):
                return 10

            def isStreamClosed(self):
                return True

            def close(self):
                self.closedCalls += 1

        class _Harness:
            finished = False
            insertedIds = set()
            processedIds = set()
            _lastInputId = 0
            isStreamClosed = False
            lastRound = False
            limitReach = False
            timerOut = False
            outputSize = _Value(10)

            def __init__(self):
                self.inputSet = _InputSet()
                self.insertedBatches = []

            def _loadInputSet(self, _):
                return self.inputSet

            def _discoverIdsAfter(self, inputSet, lastId):
                return ProtDataCounter._discoverIdsAfter(
                    self,
                    inputSet,
                    lastId,
                )

            def _reconcileClosedStreamIds(
                    self,
                    inputSet,
                    discoveredIds,
                    knownIds,
                    producerClosed,
            ):
                return ProtDataCounter._reconcileClosedStreamIds(
                    self,
                    inputSet,
                    discoveredIds,
                    knownIds,
                    producerClosed,
                )

            def isContinued(self):
                return False

            def _insertNewImageSteps(self, newIds):
                ids = list(newIds)
                self.insertedBatches.append(ids)
                self.insertedIds.update(ids)
                return []

            def info(self, message):
                pass

        protocol = _Harness()

        with patch(
                "emfacilities.protocols.protocol_data_counter.time.sleep",
                return_value=None,
        ):
            ProtDataCounter._checkNewInput(protocol)
            ProtDataCounter._checkNewInput(protocol)

        self.assertEqual(
            protocol.insertedIds,
            set(range(1, 11)),
            "Closing the stream must trigger reconciliation when the "
            "declared Set size is larger than the IDs discovered through "
            "the monotonic watermark.",
        )
        self.assertIn(
            ("id", None),
            protocol.inputSet.uniqueCalls,
            "The terminal mismatch must use a one-off full ID "
            "reconciliation instead of advancing the watermark forever.",
        )


class TestDataCounterCheckNewOutputVisibilityRegression(tests.unittest.TestCase):

    def testCheckNewOutputSkipsImageNotYetVisibleAndDoesNotFinishPrematurely(self):
        # Regression test: an id discovered earlier (via _discoverIdsAfter
        # in _checkNewInput) is not guaranteed to still be selectable via
        # Set.getItem() on a freshly reloaded input Set later - e.g. under
        # replication lag on a PostgreSQL-backed compatibility bridge.
        # Set.getItem raises rather than returning None for a missing row,
        # so a membership check is required before indexing. Skipping the
        # invisible id must not let the protocol declare itself finished
        # (and hence close the output) before that id is actually
        # persisted.
        class _Value:
            def __init__(self, value):
                self._value = value

            def get(self):
                return self._value

        class _Image:
            def __init__(self, objId):
                self.objId = objId

            def clone(self):
                return _Image(self.objId)

        class _InputSet:
            def __init__(self, visibleIds, size):
                self._visibleIds = set(visibleIds)
                self._size = size
                self.closed = False

            def getSize(self):
                return self._size

            def __contains__(self, objId):
                return objId in self._visibleIds

            def getItem(self, field, value):
                assert field == "id"
                # Real Set.getItem raises (UnboundLocalError) rather than
                # returning None for a row it cannot find - match that
                # here so a missing membership guard is caught.
                if value not in self._visibleIds:
                    raise UnboundLocalError("row not found for id %r" % value)
                return _Image(value)

            def close(self):
                self.closed = True

        class _OutputSet:
            STREAM_OPEN = 1
            STREAM_CLOSED = 2

            def __init__(self):
                self.ids = []

            def getSize(self):
                return len(self.ids)

            def append(self, image):
                self.ids.append(image.objId)

        class _Harness:
            finished = False
            processedIds = {4, 5}  # 5 is not yet visible
            isStreamClosed = True
            limitReach = False
            timerOut = False
            outputSize = _Value(100)  # far from the limit
            _inputClass = object
            _baseName = "images.sqlite"

            def __init__(self):
                # getSize() == 2 matches len(processedIds), so the
                # pre-loop optimistic completion check would (incorrectly)
                # already consider the round finished.
                self.inputSet = _InputSet(visibleIds={4}, size=2)
                self.outputSet = _OutputSet()
                self.updated = False
                self.streamMode = None
                self.errors = []

            def _loadInputSet(self, _):
                return self.inputSet

            def _loadOutputSet(self, SetClass, baseName, outputName=None):
                return self.outputSet

            def _updateOutputSet(self, outputName, outputSet, streamMode):
                self.updated = True
                self.streamMode = streamMode

            def _getFirstJoinStep(self):
                return None

            def _store(self):
                pass

            def error(self, msg):
                self.errors.append(msg)

        protocol = _Harness()

        ProtDataCounter._checkNewOutput(protocol)

        self.assertEqual([4], protocol.outputSet.ids)
        self.assertEqual({5}, protocol.processedIds)
        self.assertEqual(1, len(protocol.errors))
        self.assertFalse(
            protocol.finished,
            "A skipped (not-yet-visible) id must prevent the protocol "
            "from declaring itself finished this round.",
        )
        self.assertEqual(protocol.outputSet.STREAM_OPEN, protocol.streamMode)


class TestDataCounterPersistedOutputRegression(tests.unittest.TestCase):
    def testPersistedOutputIsRefreshedAndBackingFileIsNotWorkflowIdentity(self):
        class _Pointer:
            def __init__(self, value):
                self._value = value

            def get(self):
                return self._value

        class _ExistingOutput:
            def __init__(self):
                self.loaded = False
                self.appendEnabled = False
                self.copiedFrom = None

            def loadAllProperties(self):
                self.loaded = True

            def enableAppend(self):
                if not self.loaded:
                    raise AssertionError("Logical output must be refreshed before enableAppend().")
                self.appendEnabled = True

            def copyInfo(self, inputs):
                self.copiedFrom = inputs

            def getSize(self):
                if not self.loaded:
                    raise AssertionError("Logical output must be refreshed before reading its size.")
                return 2

            def getIdSet(self):
                if not self.loaded:
                    raise AssertionError("Logical output must be refreshed before reading its ids.")
                return {1, 2}

        class _FreshOutput:
            STREAM_OPEN = 1

            def __init__(self, filename=None):
                self.filename = filename
                self.loaded = False
                self.streamState = None
                self.copiedFrom = None

            def loadAllProperties(self):
                self.loaded = True
                raise AssertionError("A backing file must not restore an output absent from protocol outputs.")

            def setStreamState(self, state):
                self.streamState = state

            def copyInfo(self, inputs):
                self.copiedFrom = inputs

        inputs = object()
        existing = _ExistingOutput()

        class _Harness:
            def __init__(self):
                self.outputSet = existing
                self.inputImages = _Pointer(inputs)

            def _getPath(self, name):
                return "/tmp/" + name

        protocol = _Harness()

        loaded = ProtDataCounter._loadOutputSet(protocol, object, "images.sqlite", outputName=OUTPUT)
        doneIds, size = ProtDataCounter._getAllDoneIds(protocol)

        self.assertIs(existing, loaded)
        self.assertTrue(existing.loaded)
        self.assertTrue(existing.appendEnabled)
        self.assertEqual({1, 2}, set(doneIds))
        self.assertEqual(2, size)

        del protocol.outputSet

        with patch("emfacilities.protocols.protocol_data_counter.pwutils.cleanPath") as cleanPathMock:
            fresh = ProtDataCounter._loadOutputSet(protocol, _FreshOutput, "images.sqlite", outputName=OUTPUT)

        cleanPathMock.assert_called_once_with("/tmp/images.sqlite")
        self.assertFalse(fresh.loaded)
        self.assertEqual(_FreshOutput.STREAM_OPEN, fresh.streamState)
        self.assertIs(inputs, fresh.copiedFrom)
