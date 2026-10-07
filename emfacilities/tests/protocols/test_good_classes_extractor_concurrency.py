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
"""extractElements runs under STEPS_PARALLEL, on several threads at once.

Everything it touches that is shared - the output Sets, the particle
counters, and the pyplot global state the plots are drawn on - has to be
reached under the protocol lock. pyplot in particular keeps a single
"current figure" per process, so two steps drawing at the same time can
save each other's figure.
"""
import threading
import time
import unittest

from emfacilities.protocols.protocol_good_classes_extractor import (
    ProtGoodClassesExtractor,
    OUTPUT_PARTICLES,
    OUTPUT_DISCARDED_PARTICLES,
)


class _TrackingLock:
    """A lock that remembers how many threads are inside it."""

    def __init__(self):
        self._lock = threading.Lock()
        self.heldBy = set()

    def __enter__(self):
        self._lock.acquire()
        self.heldBy.add(threading.current_thread().ident)
        return self

    def __exit__(self, *exc):
        self.heldBy.discard(threading.current_thread().ident)
        self._lock.release()
        return False

    def heldHere(self):
        return threading.current_thread().ident in self.heldBy


class _Image:
    def __init__(self, objId, creation):
        self._objId = objId
        self._creation = creation

    def getObjId(self):
        return self._objId

    def getObjCreation(self):
        return self._creation

    def clone(self):
        return _Image(self._objId, self._creation)


class _Class:
    def __init__(self, objId, images):
        self._objId = objId
        self._images = images

    def getObjId(self):
        return self._objId

    def iterItems(self, orderBy=None, direction=None, where=None):
        return iter(self._images)


class _ClassesStub:
    def __init__(self, classes):
        self._classes = classes

    def iterItems(self, orderBy=None, direction=None, where=None):
        return iter(self._classes)


class _OutputStub:
    def __init__(self):
        self.items = []

    def append(self, item):
        self.items.append(item)

    def getIdSet(self):
        return {i.getObjId() for i in self.items}

    def __len__(self):
        return len(self.items)


class _ExtractHarness:
    """Stand-in exercising the real extractElements body."""

    def __init__(self):
        self._lock = _TrackingLock()
        self.dictsTimes = {}
        self.goodClassesIDs = [1]
        self.goodParticles = []
        self.badParticles = []
        self._processedParticleIds = set()
        self.particlesDistribution = {'good': [], 'bad': []}
        self.isStreamClosed = False
        self.plotCalls = 0
        self.plotsOutsideTheLock = []
        self.concurrentPlotters = 0
        self.maxConcurrentPlotters = 0
        self._counterLock = threading.Lock()
        self._outputs = {
            OUTPUT_PARTICLES: _OutputStub(),
            OUTPUT_DISCARDED_PARTICLES: _OutputStub(),
        }

    def _loadOutputSet(self, outputName, suffix):
        return self._outputs[outputName]

    def _updateOutputSet(self, outputName, outputSet, state):
        setattr(self, outputName, outputSet)

    def _createPlots(self):
        # Record how this was reached, then hold the slot long enough for
        # a second thread to collide with it.
        with self._counterLock:
            self.plotCalls += 1
            self.plotsOutsideTheLock.append(not self._lock.heldHere())
            self.concurrentPlotters += 1
            self.maxConcurrentPlotters = max(
                self.maxConcurrentPlotters, self.concurrentPlotters
            )

        # Hold the slot open long enough that an unsynchronised second
        # thread is bound to walk into it.
        time.sleep(0.05)

        with self._counterLock:
            self.concurrentPlotters -= 1

    def info(self, *args):
        pass

    def debug(self, *args):
        pass


class TestExtractElementsPlotsUnderTheLock(unittest.TestCase):

    def _buildHarness(self):
        return _ExtractHarness()

    def _run(self, harness, classes):
        ProtGoodClassesExtractor.extractElements(
            harness, _ClassesStub(classes)
        )

    def testPlottingHappensUnderTheProtocolLock(self):
        harness = self._buildHarness()

        self._run(harness, [_Class(1, [_Image(10, 'c1')])])

        self.assertEqual(harness.plotCalls, 1)
        self.assertEqual(
            harness.plotsOutsideTheLock,
            [False],
            "The plots are drawn on pyplot's single global current "
            "figure: reaching them without the protocol lock lets two "
            "parallel steps save each other's figure.",
        )

    def testTwoParallelStepsNeverPlotAtTheSameTime(self):
        harness = self._buildHarness()
        errors = []

        def runStep(classId, imageId):
            try:
                self._run(
                    harness, [_Class(classId, [_Image(imageId, 'c1')])]
                )
            except Exception as exc:  # noqa: BLE001 - reported below
                errors.append(exc)

        threads = [
            threading.Thread(target=runStep, args=(1, 10)),
            threading.Thread(target=runStep, args=(1, 11)),
        ]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=10)

        for thread in threads:
            self.assertFalse(
                thread.is_alive(), "extractElements deadlocked: %r" % errors
            )

        self.assertEqual(
            harness.maxConcurrentPlotters,
            1,
            "Two parallel steps drew plots at the same time: pyplot keeps "
            "one global current figure per process, so one step's savefig "
            "can write out the other step's half-drawn figure.",
        )


if __name__ == '__main__':
    unittest.main()
