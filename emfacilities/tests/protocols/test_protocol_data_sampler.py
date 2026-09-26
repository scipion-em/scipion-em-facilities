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
    def getIdSet(self):
        return {1, 2, 3, 4, 5, 6}

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
            def getIdSet(self):
                return {1, 2, 3, 4, 5, 6}

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
        self.assertEqual(prot.sampleIds, {2})


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
        prot.sampleIds = {1, 4}
        prot.inputFn = 'input.sqlite'
        prot._inputClass = object
        prot._baseName = 'images.sqlite'
        prot._getAllDoneIds = lambda: ([1, 4], 2)
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
                self.enableAppendCalls = 0
                self.copiedFrom = None

            def enableAppend(self):
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
        self.assertEqual(1, existingOutputSet.enableAppendCalls)

class TestDataSamplerFinalizationRegression(tests.unittest.TestCase):
    def testFinishedStepsCheckDoesNotTouchInputOrOutput(self):
        class _Harness:
            finished = True

            def __init__(self):
                self._checkNewInput = Mock()
                self._checkNewOutput = Mock()

        protocol = _Harness()

        ProtDataSampler._stepsCheck(protocol)

        protocol._checkNewInput.assert_not_called()
        protocol._checkNewOutput.assert_not_called()
