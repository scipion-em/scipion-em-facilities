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
"""Resuming must see the steps the previous run actually finished.

pyworkflow carries a previous run's steps in _prevSteps and only copies
the ones whose index also exists in the freshly inserted _steps. A steps
generator inserts its work while it runs, long after that comparison, so
on Resume _steps holds the generator alone and every finished sampling
step of the previous run exists only in _prevSteps. Reading _steps by
itself therefore restores nothing at all.
"""
import json
import os
import shutil
import tempfile
import unittest

from emfacilities.protocols.protocol_data_sampler import ProtDataSampler
from emfacilities.protocols.protocol_streaming_base import (
    ProtFacilitiesStreamingBase,
)


class _Value:
    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value


class _Step:
    def __init__(self, funcName, ids=None, finished=True, stateFile=None):
        self.funcName = _Value(funcName)
        self._finished = finished
        self.argsStr = _Value(
            json.dumps([ids]) if ids is not None else None
        )
        self._resultFiles = _Value(
            json.dumps([stateFile]) if stateFile else None
        )

    def isFinished(self):
        return self._finished


class _Harness:
    _restoreRuntimeStateFromFinishedSteps = (
        ProtDataSampler._restoreRuntimeStateFromFinishedSteps
    )
    _iterKnownSteps = ProtFacilitiesStreamingBase._iterKnownSteps

    def __init__(self, steps=(), prevSteps=()):
        self.insertedIds = set()
        self.processedIds = set()
        self.sampleIds = set()
        self._steps = list(steps)
        self._prevSteps = list(prevSteps)


class TestResumeReadsPreviousRunSteps(unittest.TestCase):

    def setUp(self):
        self.workdir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.workdir, True)

    def _stateFile(self, name, processedIds, sampledIds):
        path = os.path.join(self.workdir, name)
        with open(path, 'w', encoding='utf-8') as handle:
            json.dump(
                {"processedIds": processedIds, "sampledIds": sampledIds},
                handle,
            )
        return path

    def testFinishedStepsOnlyInPrevStepsAreRestored(self):
        stateFile = self._stateFile("batch.json", [1, 2, 3, 4], [2, 4])
        harness = _Harness(
            # What a resumed run starts with: the generator, nothing else.
            steps=[_Step("stepsGeneratorStep", finished=False)],
            prevSteps=[
                _Step("stepsGeneratorStep", finished=False),
                _Step("samplingStep", ids=[1, 2, 3, 4],
                      stateFile=stateFile),
            ],
        )

        harness._restoreRuntimeStateFromFinishedSteps(doneIds=set())

        self.assertEqual(
            harness.processedIds,
            {1, 2, 3, 4},
            "The previous run sampled this batch: not seeing it makes the "
            "resumed run sample the very same images all over again.",
        )
        self.assertEqual(
            harness.sampleIds,
            {2, 4},
            "These images were selected but never written out. The "
            "selection is random, so nothing can recover them except the "
            "record the finished step left behind.",
        )

    def testStepsAlreadyCopiedIntoStepsAreNotCountedTwice(self):
        stateFile = self._stateFile("batch.json", [1, 2], [2])
        step = _Step("samplingStep", ids=[1, 2], stateFile=stateFile)
        harness = _Harness(steps=[step], prevSteps=[step])

        harness._restoreRuntimeStateFromFinishedSteps(doneIds=set())

        self.assertEqual(harness.processedIds, {1, 2})
        self.assertEqual(harness.sampleIds, {2})

    def testAlreadyPersistedSelectionsAreNotQueuedAgain(self):
        stateFile = self._stateFile("batch.json", [1, 2, 3, 4], [2, 4])
        harness = _Harness(
            steps=[_Step("stepsGeneratorStep", finished=False)],
            prevSteps=[
                _Step("samplingStep", ids=[1, 2, 3, 4],
                      stateFile=stateFile),
            ],
        )

        # Image 2 already reached the output before the run stopped.
        harness._restoreRuntimeStateFromFinishedSteps(doneIds={2})

        self.assertEqual(
            harness.sampleIds,
            {4},
            "Image 2 is durable output already: queueing it again would "
            "append it to the output a second time.",
        )

    def testUnfinishedStepsAreIgnored(self):
        stateFile = self._stateFile("batch.json", [1, 2], [2])
        harness = _Harness(
            steps=[],
            prevSteps=[
                _Step("samplingStep", ids=[1, 2], stateFile=stateFile,
                      finished=False),
            ],
        )

        harness._restoreRuntimeStateFromFinishedSteps(doneIds=set())

        self.assertEqual(harness.processedIds, set())
        self.assertEqual(harness.sampleIds, set())

    def testAMissingPrevStepsAttributeIsTolerated(self):
        harness = _Harness(steps=[])
        del harness._prevSteps

        harness._restoreRuntimeStateFromFinishedSteps(doneIds={7})

        self.assertEqual(harness.processedIds, {7})


if __name__ == '__main__':
    unittest.main()
