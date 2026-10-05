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
import os
import json
from unittest.mock import Mock, patch
from pyworkflow.tests import BaseTest, DataSet
from pwem.protocols.protocol_import import ProtImportMicrographs
from pyworkflow.object import Pointer
import pyworkflow.tests as tests
from emfacilities.protocols.protocol_data_sampler import ProtDataSampler, OUTPUT



class _SamplerValue:
    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value

    def hasValue(self):
        return self.value not in (None, "")

    def __eq__(self, other):
        return self.value == other


class _SamplerInputSet:
    def getUniqueValues(self, attributes, where=None):
        # Continue starts discovery again from watermark 0 and reconciles
        # these discovered IDs with persisted/runtime state.
        return [1, 2, 3, 4, 5, 6]

    def isStreamClosed(self):
        return False

    def close(self):
        pass


class _FinishedSamplingStep:
    funcName = _SamplerValue("samplingStep")
    argsStr = _SamplerValue("[[1, 2, 3]]")

    def __init__(self, resultFile=None):
        if resultFile is not None:
            self._resultFiles = _SamplerValue(
                json.dumps([resultFile])
            )

    def isFinished(self):
        return True


class TestDataSampler(BaseTest):

    def _prepareContinueProtocol(self, steps, doneIds):
        insertedBatches = []

        prot = self.newProtocol(
            ProtDataSampler,
            batchSize=3,
            samplingProportion=0.5,
        )
        prot.inputFn = "input.sqlite"
        prot.insertedIds = set()
        prot.processedIds = set()
        prot.sampleIds = set()
        prot._lastInputId = 0
        prot._pendingInputIds = []
        prot.isStreamClosed = False
        prot._steps = steps
        prot.isContinued = lambda: True
        prot._loadInputSet = lambda _: _SamplerInputSet()
        prot._getAllDoneIds = lambda: (list(doneIds), len(doneIds))
        prot._getFirstJoinStep = lambda: None
        prot.updateSteps = lambda: None

        def insertNewImageSteps(newIds, batchSize):
            insertedBatches.append(list(newIds))
            return []

        prot._insertNewImageSteps = insertNewImageSteps
        return prot, insertedBatches

    """ Test data sampler protocol """

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


    def testDataSampler25(self):
        prot = self._runDataSampler("Random sampling of 0.35 proportion", batch=4, proportion=0.35)
        self.assertSetSize(prot.outputSet, size=1)


    def testDataSampler50(self):
        prot = self._runDataSampler("Random sampling of 0.70 proportion", batch=4,  proportion=0.70)
        self.assertSetSize(prot.outputSet, size=2)




    def testSamplingStepPersistsChosenIdsForContinue(self):
        prot = self.newProtocol(
            ProtDataSampler,
            batchSize=3,
            samplingProportion=0.5,
        )
        prot.processedIds = set()
        prot.sampleIds = set()

        with patch(
                "emfacilities.protocols.protocol_data_sampler.sample_proportion",
                return_value=[2],
        ):
            stateFile = prot.samplingStep([1, 2, 3])

        self.assertTrue(stateFile)
        self.assertTrue(os.path.exists(stateFile))

        with open(stateFile, "r", encoding="utf-8") as handle:
            state = json.load(handle)

        self.assertEqual(state["processedIds"], [1, 2, 3])
        self.assertEqual(state["sampledIds"], [2])
        self.assertEqual(prot.processedIds, {1, 2, 3})
        self.assertEqual(prot.sampleIds, {2})



    def testDetectsNewInputWhenMtimeDoesNotChange(self):
        class InputSet:
            def getUniqueValues(self, attributes, where=None):
                self.query = (attributes, where)
                return [4, 5, 6]

            def isStreamClosed(self):
                return False

            def close(self):
                pass

        insertedBatches = []

        prot = self.newProtocol(
            ProtDataSampler,
            batchSize=3,
            samplingProportion=0.5,
        )
        prot.inputFn = "input.sqlite"
        prot.insertedIds = {1, 2, 3}
        prot.processedIds = {1, 2, 3}
        prot.sampleIds = {2}
        prot._lastInputId = 3
        prot._pendingInputIds = []
        prot.isStreamClosed = False
        prot.isContinued = lambda: False
        prot._loadInputSet = lambda _: InputSet()
        prot._getFirstJoinStep = lambda: None
        prot.updateSteps = lambda: None

        def insertNewImageSteps(newIds, batchSize):
            insertedBatches.append(list(newIds))
            return []

        prot._insertNewImageSteps = insertNewImageSteps

        with patch(
                "emfacilities.protocols.protocol_data_sampler.os.path.getmtime",
                return_value=0,
        ):
            prot._checkNewInput()

        self.assertEqual(
            insertedBatches,
            [[4, 5, 6]],
            "New logical Set items must be detected even when the sqlite "
            "mtime does not change.",
        )


    def testContinueRestoresSampleChosenBeforeOutputFlush(self):
        prot = self.newProtocol(
            ProtDataSampler,
            batchSize=3,
            samplingProportion=0.5,
        )
        stateFile = prot._getExtraPath("sampling_state_test.json")
        os.makedirs(os.path.dirname(stateFile), exist_ok=True)

        with open(stateFile, "w", encoding="utf-8") as handle:
            json.dump(
                {
                    "processedIds": [1, 2, 3],
                    "sampledIds": [2],
                },
                handle,
            )

        prot, insertedBatches = self._prepareContinueProtocol(
            [_FinishedSamplingStep(stateFile)],
            doneIds=[],
        )

        with patch(
                "emfacilities.protocols.protocol_data_sampler.os.path.getmtime",
                return_value=0,
        ):
            prot._checkNewInput()

        self.assertEqual(insertedBatches, [[4, 5, 6]])
        self.assertEqual(prot.insertedIds, {1, 2, 3})
        self.assertEqual(prot.processedIds, {1, 2, 3})
        self.assertEqual(prot.sampleIds, {2})


    def testContinueRestoresFinishedSamplingBatches(self):
        prot, insertedBatches = self._prepareContinueProtocol(
            [_FinishedSamplingStep()],
            doneIds=[2],
        )

        with patch(
                "emfacilities.protocols.protocol_data_sampler.os.path.getmtime",
                return_value=0,
        ):
            prot._checkNewInput()

        self.assertEqual(insertedBatches, [[4, 5, 6]])
        self.assertEqual(prot.insertedIds, {1, 2, 3})
        self.assertEqual(prot.processedIds, {1, 2, 3})
        self.assertEqual(
            prot.sampleIds,
            set(),
            "Persisted doneIds must not be restored as pending samples.",
        )


    def testClosedStreamFinishesAfterAllInputWasProcessed(self):
        class InputSet:
            def getSize(self):
                return 6

            def getIdSet(self):
                return {1, 2, 3, 4, 5, 6}

            def close(self):
                pass

        prot = self.newProtocol(ProtDataSampler, batchSize=3, samplingProportion=0.5)
        prot.finished = False
        prot.isStreamClosed = True
        prot.processedIds = {1, 2, 3, 4, 5, 6}
        # IDs 1 and 4 are already persisted output; sampleIds contains only
        # samples still pending persistence.
        prot.sampleIds = set()
        prot.inputFn = 'input.sqlite'
        prot._inputClass = object
        prot._baseName = 'images.sqlite'
        prot._loadInputSet = lambda _: InputSet()
        prot._loadOutputSet = lambda *args, **kwargs: object()
        prot._updateOutputSet = lambda *args, **kwargs: None
        prot._getFirstJoinStep = lambda: None
        prot._store = lambda: None

        prot._checkNewOutput()

        self.assertTrue(prot.finished)




    def _runDataSampler(cls, label, batch, proportion):
        protDataSampler = cls.newProtocol(ProtDataSampler,
                                          batchSize=batch,
                                          samplingProportion=proportion,
                                          delay=3)
        protDataSampler.inputImages = Pointer(cls.protImport, extended='outputMicrographs')
        protDataSampler.setObjLabel(label)
        cls.launchProtocol(protDataSampler)

        return protDataSampler


    def testSamplingStepPublishesCompletionAfterSampleState(self):
        prot = self.newProtocol(
            ProtDataSampler,
            batchSize=3,
            samplingProportion=0.5,
        )
        prot.sampleIds = set()

        stateFile = prot._getSamplingStateFile([1, 2, 3])
        if os.path.exists(stateFile):
            os.remove(stateFile)

        class _ProcessedIds(set):
            def update(innerSelf, values):
                self.assertEqual(
                    {2},
                    prot.sampleIds,
                    "sampleIds must be published before processedIds marks "
                    "the batch complete.",
                )
                self.assertTrue(
                    os.path.exists(stateFile),
                    "The sampling state must be durable before processedIds "
                    "marks the batch complete.",
                )
                super(_ProcessedIds, innerSelf).update(values)

        prot.processedIds = _ProcessedIds()

        with patch(
                "emfacilities.protocols.protocol_data_sampler.sample_proportion",
                return_value=[2],
        ):
            prot.samplingStep([1, 2, 3])

        self.assertEqual({1, 2, 3}, prot.processedIds)
        self.assertEqual({2}, prot.sampleIds)


class TestDataSamplerLoadOutputSet(tests.unittest.TestCase):
    """Lightweight regression tests that need no real project/dataset."""

    def testLoadOutputSetReusesLogicalOutputWithoutBackingFile(self):
        # Regression test: an output that Scipion already knows about
        # (protocol.outputSet) must be reused even when its backing file
        # was never materialized on disk yet. Falling through to "no
        # backing file -> build a fresh, empty Set" would silently discard
        # whatever was already appended to the real logical output.
        prot = ProtDataSampler()

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
                "emfacilities.protocols.protocol_data_sampler.os.path.exists",
                return_value=False,
        ):
            outputSet = prot._loadOutputSet(
                object, "images.sqlite", outputName=OUTPUT
            )

        self.assertIs(existingOutputSet, outputSet)
        self.assertEqual(1, existingOutputSet.loadAllPropertiesCalls)
        self.assertEqual(1, existingOutputSet.enableAppendCalls)

class TestDataSamplerFinalizationRegression(tests.unittest.TestCase):
    def testFinishedGeneratorDoesNotPollInputOrOutput(self):
        class _Harness:
            def __init__(self):
                self.finished = False
                self._checkNewInput = Mock()
                self._checkNewOutput = Mock()
                self._closeOutputSet = Mock()

            def initializeParams(self):
                self.finished = True

        protocol = _Harness()

        ProtDataSampler.stepsGeneratorStep(protocol)

        protocol._checkNewInput.assert_not_called()
        protocol._checkNewOutput.assert_not_called()
        protocol._closeOutputSet.assert_called_once_with()


class TestDataSamplerBackendIndependence(tests.unittest.TestCase):
    def testInputSetAccessUsesLogicalPointerWithoutBackingFilename(self):
        class _ForbiddenSetConstructor:
            def __init__(self, *args, **kwargs):
                raise AssertionError(
                    "DataSampler must not reconstruct the input Set from a backing filename."
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
                    "DataSampler must not require a backing filename for the logical input Set."
                )

            def loadAllProperties(self):
                self.loadCalls += 1

        class _Pointer:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        logicalInput = _LogicalInputSet()
        protocol = ProtDataSampler()
        protocol.inputImages = _Pointer(logicalInput)

        protocol.initializeParams()

        loadedInput = protocol._loadInputSet("must-not-be-used.sqlite")

        self.assertIs(logicalInput, loadedInput)
        self.assertEqual(1, logicalInput.loadCalls)


class TestDataSamplerStreamingScalability(tests.unittest.TestCase):
    def testIncrementalDiscoveryKeepsPartialBatchAcrossPolls(self):
        class _Value:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        class _InputSet:
            def __init__(self):
                self.poll = 0
                self.uniqueCalls = []
                self.closed = 0

            def getIdSet(self):
                raise AssertionError(
                    "DataSampler must not fetch the complete input ID set on every poll."
                )

            def getUniqueValues(self, attributes, where=None):
                self.uniqueCalls.append((attributes, where))
                self.poll += 1
                if self.poll == 1:
                    return [4, 5]
                if self.poll == 2:
                    return [6]
                return []

            def isStreamClosed(self):
                return False

            def close(self):
                self.closed += 1

        class _Harness:
            insertedIds = {1, 2, 3}
            processedIds = {1, 2, 3}
            sampleIds = {2}
            isStreamClosed = False
            _lastInputId = 3
            _pendingInputIds = []
            batchSize = _Value(3)

            def __init__(self):
                self.inputSet = _InputSet()
                self.insertedBatches = []

            def _loadInputSet(self, _):
                return self.inputSet

            def _discoverIdsAfter(self, inputSet, lastId):
                return ProtDataSampler._discoverIdsAfter(
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
                return ProtDataSampler._reconcileClosedStreamIds(
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

            def _insertNewImageSteps(self, newIds, batchSize):
                self.insertedBatches.append(list(newIds))
                self.insertedIds.update(newIds)
                return []

            def updateSteps(self):
                pass

        protocol = _Harness()

        ProtDataSampler._checkNewInput(protocol)

        self.assertEqual([], protocol.insertedBatches)
        self.assertEqual([4, 5], protocol._pendingInputIds)
        self.assertEqual(5, protocol._lastInputId)

        ProtDataSampler._checkNewInput(protocol)

        self.assertEqual([[4, 5, 6]], protocol.insertedBatches)
        self.assertEqual([], protocol._pendingInputIds)
        self.assertEqual(6, protocol._lastInputId)
        self.assertEqual(
            [("id", "id > 3"), ("id", "id > 5")],
            protocol.inputSet.uniqueCalls,
        )
        self.assertEqual(2, protocol.inputSet.closed)


class TestDataSamplerOutputPollingScalability(tests.unittest.TestCase):
    def testOutputPollingUsesSizesAndPendingSamples(self):
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

            def getIdSet(self):
                raise AssertionError(
                    "_checkNewOutput must not materialize every input ID."
                )

            def getItem(self, field, value):
                return _Image(value)

            def close(self):
                self.closed = True

        class _OutputSet:
            STREAM_OPEN = 1
            STREAM_CLOSED = 2

            def __init__(self):
                self.ids = [2]

            def getSize(self):
                return len(self.ids)

            def getIdSet(self):
                raise AssertionError(
                    "_checkNewOutput must not rescan every persisted output ID."
                )

            def append(self, image):
                self.ids.append(image.objId)

        class _Harness:
            finished = False
            isStreamClosed = True
            processedIds = {1, 2, 3, 4, 5}
            # Pending sampled IDs only: ID 2 is already persisted.
            sampleIds = {4, 5}
            _inputClass = object
            _baseName = "images.sqlite"

            def __init__(self):
                self.inputSet = _InputSet()
                self.outputSet = _OutputSet()
                self.updatedMode = None

            def _getAllDoneIds(self):
                return ProtDataSampler._getAllDoneIds(self)

            def _loadInputSet(self, _):
                return self.inputSet

            def _loadOutputSet(self, SetClass, baseName, outputName=None):
                return self.outputSet

            def _updateOutputSet(self, outputName, outputSet, streamMode):
                self.updatedMode = streamMode

            def _getFirstJoinStep(self):
                return None

            def _store(self):
                pass

        protocol = _Harness()

        ProtDataSampler._checkNewOutput(protocol)

        self.assertTrue(protocol.finished)
        self.assertEqual([2, 4, 5], protocol.outputSet.ids)
        self.assertEqual(set(), protocol.sampleIds)
        self.assertEqual(protocol.outputSet.STREAM_CLOSED, protocol.updatedMode)
        self.assertTrue(protocol.inputSet.closed)


class TestDataSamplerStreamingArchitecture(tests.unittest.TestCase):
    def testUsesFacilitiesStreamingBaseAndCoreGeneratorInsertion(self):
        from pyworkflow.protocol import ProtStreamingBase, STEPS_PARALLEL
        from emfacilities.protocols.protocol_streaming_base import (
            ProtFacilitiesStreamingBase,
        )

        self.assertTrue(
            issubclass(ProtDataSampler, ProtFacilitiesStreamingBase)
        )
        self.assertTrue(
            issubclass(ProtFacilitiesStreamingBase, ProtStreamingBase)
        )
        self.assertIs(
            ProtDataSampler._insertAllSteps,
            ProtStreamingBase._insertAllSteps,
        )
        self.assertEqual(
            STEPS_PARALLEL,
            ProtDataSampler.stepsExecutionMode,
        )


    def testGeneratorWaitsForPersistedCompletionAfterInputCloses(self):
        class _Harness:
            def __init__(self):
                self.finished = False
                self.isStreamClosed = False
                self.iteration = 0
                self.events = []

            def initializeParams(self):
                self.finished = False
                self.events.append("initialize")

            def _checkNewInput(self):
                self.iteration += 1
                self.isStreamClosed = True
                self.events.append("input-%d-closed" % self.iteration)

            def _checkNewOutput(self):
                self.events.append("output-%d" % self.iteration)

                # Input being closed is not enough. Simulate one extra round
                # before output persistence/completion is confirmed.
                if self.iteration == 2:
                    self.finished = True

            def _streamingSleepOnWait(self):
                self.events.append("sleep-%d" % self.iteration)

            def _closeOutputSet(self):
                self.events.append("close")

        protocol = _Harness()

        ProtDataSampler.stepsGeneratorStep(protocol)

        self.assertEqual(
            [
                "initialize",
                "input-1-closed",
                "output-1",
                "sleep-1",
                "input-2-closed",
                "output-2",
                "close",
            ],
            protocol.events,
        )

class TestDataSamplerClosedStreamLateVisibilityRegression(
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
            insertedIds = set()
            processedIds = set()
            sampleIds = set()
            _lastInputId = 0
            _pendingInputIds = []
            isStreamClosed = False
            batchSize = _Value(10)

            def __init__(self):
                self.inputSet = _InputSet()
                self.insertedBatches = []

            def _loadInputSet(self, _):
                return self.inputSet

            def _discoverIdsAfter(self, inputSet, lastId):
                return ProtDataSampler._discoverIdsAfter(
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
                return ProtDataSampler._reconcileClosedStreamIds(
                    self,
                    inputSet,
                    discoveredIds,
                    knownIds,
                    producerClosed,
                )

            def isContinued(self):
                return False

            def _insertNewImageSteps(self, newIds, batchSize):
                ids = list(newIds)
                self.insertedBatches.append(ids)
                self.insertedIds.update(ids)
                return []

            def info(self, message):
                pass

        protocol = _Harness()

        ProtDataSampler._checkNewInput(protocol)
        ProtDataSampler._checkNewInput(protocol)

        self.assertEqual(
            protocol.insertedIds,
            set(range(1, 11)),
            "A closed sampler input must reconcile IDs that became visible "
            "below an already advanced watermark.",
        )
        self.assertIn(
            ("id", None),
            protocol.inputSet.uniqueCalls,
            "The terminal mismatch must use a full ID reconciliation.",
        )


class TestDataSamplerCheckNewOutputVisibilityRegression(tests.unittest.TestCase):

    def testCheckNewOutputSkipsSampleNotYetVisibleAndDoesNotFinishPrematurely(self):
        # Regression test: a sampled id is not guaranteed to still be
        # selectable via Set.getItem() on a freshly reloaded input Set
        # later - e.g. under replication lag on a PostgreSQL-backed
        # compatibility bridge. Set.getItem raises rather than returning
        # None for a missing row, so a membership check is required
        # before indexing. self.finished here is driven only by
        # processedIds/inputSize (unrelated to whether the sample was
        # actually persisted), so skipping the invisible id must force
        # self.finished back to False - otherwise the outer generator
        # loop would exit and close the output before that sample is
        # ever persisted.
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
            isStreamClosed = True
            processedIds = {1, 2, 3, 4, 5}  # all scanned -> optimistically done
            sampleIds = {4, 5}  # 5 is not yet visible
            _inputClass = object
            _baseName = "images.sqlite"

            def __init__(self):
                self.inputSet = _InputSet(visibleIds={4}, size=5)
                self.outputSet = _OutputSet()
                self.updatedMode = None
                self.errors = []

            def _loadInputSet(self, _):
                return self.inputSet

            def _loadOutputSet(self, SetClass, baseName, outputName=None):
                return self.outputSet

            def _updateOutputSet(self, outputName, outputSet, streamMode):
                self.updatedMode = streamMode

            def _getFirstJoinStep(self):
                return None

            def _store(self):
                pass

            def error(self, msg):
                self.errors.append(msg)

        protocol = _Harness()

        ProtDataSampler._checkNewOutput(protocol)

        self.assertEqual([4], protocol.outputSet.ids)
        self.assertEqual({5}, protocol.sampleIds)
        self.assertEqual(1, len(protocol.errors))
        self.assertFalse(
            protocol.finished,
            "A skipped (not-yet-visible) sample must prevent the "
            "protocol from declaring itself finished this round.",
        )
        self.assertEqual(protocol.outputSet.STREAM_OPEN, protocol.updatedMode)


class TestDataSamplerPersistedOutputRegression(tests.unittest.TestCase):
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

        loaded = ProtDataSampler._loadOutputSet(protocol, object, "images.sqlite", outputName=OUTPUT)
        doneIds, size = ProtDataSampler._getAllDoneIds(protocol)

        self.assertIs(existing, loaded)
        self.assertTrue(existing.loaded)
        self.assertTrue(existing.appendEnabled)
        self.assertEqual({1, 2}, set(doneIds))
        self.assertEqual(2, size)

        del protocol.outputSet

        with patch("emfacilities.protocols.protocol_data_sampler.pwutils.cleanPath") as cleanPathMock:
            fresh = ProtDataSampler._loadOutputSet(protocol, _FreshOutput, "images.sqlite", outputName=OUTPUT)

        cleanPathMock.assert_called_once_with("/tmp/images.sqlite")
        self.assertFalse(fresh.loaded)
        self.assertEqual(_FreshOutput.STREAM_OPEN, fresh.streamState)
        self.assertIs(inputs, fresh.copiedFrom)
