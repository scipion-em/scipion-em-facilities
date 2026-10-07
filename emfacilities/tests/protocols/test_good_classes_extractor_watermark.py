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
"""The extractor advances on creation time, which is stored by the second.

pyworkflow stores an object's creation stamp in UTC with no microseconds
(see SqliteFlatMapper.fmtDate), so a batch of particles written within
the same second all carry the exact same value. A watermark that only
accepts strictly later stamps therefore drops every sibling of the last
particle it happened to see - silently, and for good.
"""
import threading
import unittest

from emfacilities.protocols.protocol_good_classes_extractor import (
    ProtGoodClassesExtractor,
    OUTPUT_PARTICLES,
    OUTPUT_DISCARDED_PARTICLES,
)


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
    """A class whose visible particles can grow between polls."""

    def __init__(self, objId, images):
        self._objId = objId
        self.images = list(images)
        self.queries = []

    def getObjId(self):
        return self._objId

    def iterItems(self, orderBy=None, direction=None, where=None):
        self.queries.append(where)
        selected = sorted(self.images, key=lambda i: (i.getObjCreation(),
                                                      i.getObjId()))
        if where is None:
            return iter(selected)

        # Mimic the comparison the Set would run server side.
        if where.startswith('creation>='):
            bound = where.split('"')[1]
            return iter([i for i in selected
                         if i.getObjCreation() >= bound])
        if where.startswith('creation>'):
            bound = where.split('"')[1]
            return iter([i for i in selected
                         if i.getObjCreation() > bound])

        raise AssertionError("Unexpected query: %r" % (where,))


class _ClassesStub:
    def __init__(self, classes):
        self.classes = classes
        self.closed = 0

    def iterItems(self, orderBy=None, direction=None, where=None):
        return iter(self.classes)

    def getStreamState(self):
        return False

    def close(self):
        self.closed += 1


class _OutputStub:
    def __init__(self):
        self.items = []
        self.idSetCalls = 0

    def append(self, item):
        self.items.append(item)

    def getIdSet(self):
        self.idSetCalls += 1
        return {i.getObjId() for i in self.items}

    def __len__(self):
        return len(self.items)


class _Harness:
    """Stand-in exercising the real extractor bodies."""

    _newParticlesToProcess = ProtGoodClassesExtractor._newParticlesToProcess

    def __init__(self, classes):
        self._lock = threading.Lock()
        self.classesStub = _ClassesStub(classes)
        self.dictsTimes = {}
        self.goodClassesIDs = [1]
        self.goodParticles = []
        self.badParticles = []
        self._processedParticleIds = set()
        self.particlesDistribution = {'good': [], 'bad': []}
        self.isStreamClosed = False
        self._outputs = {
            OUTPUT_PARTICLES: _OutputStub(),
            OUTPUT_DISCARDED_PARTICLES: _OutputStub(),
        }

    def _loadInputClassesSet(self):
        return self.classesStub

    def _loadOutputSet(self, outputName, suffix):
        return self._outputs[outputName]

    def _updateOutputSet(self, outputName, outputSet, state):
        setattr(self, outputName, outputSet)

    def _createPlots(self):
        pass

    def info(self, *args):
        pass

    def debug(self, *args):
        pass

    def extract(self):
        ProtGoodClassesExtractor.extractElements(self, self.classesStub)


class TestCreationWatermarkKeepsSameSecondSiblings(unittest.TestCase):

    def testParticlesSharingTheWatermarkSecondAreNotLost(self):
        # Three particles written within the same second; only the first
        # two are visible when the first extraction runs.
        clazz = _Class(1, [_Image(10, '2026-01-01 10:00:00'),
                           _Image(11, '2026-01-01 10:00:00')])
        harness = _Harness([clazz])

        harness.extract()

        clazz.images.append(_Image(12, '2026-01-01 10:00:00'))
        harness.extract()

        collected = sorted(
            i.getObjId() for i in harness._outputs[OUTPUT_PARTICLES].items
        )
        self.assertEqual(
            collected,
            [10, 11, 12],
            "Particle 12 shares its creation second with the watermark: a "
            "strictly-greater comparison drops it for good, and creation "
            "stamps are stored without microseconds.",
        )

    def testAlreadyCollectedParticlesAreNotDuplicated(self):
        clazz = _Class(1, [_Image(10, '2026-01-01 10:00:00'),
                           _Image(11, '2026-01-01 10:00:00')])
        harness = _Harness([clazz])

        harness.extract()
        harness.extract()
        harness.extract()

        collected = sorted(
            i.getObjId() for i in harness._outputs[OUTPUT_PARTICLES].items
        )
        self.assertEqual(
            collected,
            [10, 11],
            "Re-reading the boundary second must not append the same "
            "particle again.",
        )

    def testDiscardedParticlesKeepTheSameGuarantee(self):
        clazz = _Class(7, [_Image(10, '2026-01-01 10:00:00'),
                           _Image(11, '2026-01-01 10:00:00')])
        harness = _Harness([clazz])
        harness.goodClassesIDs = [1]  # class 7 is discarded

        harness.extract()
        clazz.images.append(_Image(12, '2026-01-01 10:00:00'))
        harness.extract()

        collected = sorted(
            i.getObjId()
            for i in harness._outputs[OUTPUT_DISCARDED_PARTICLES].items
        )
        self.assertEqual(collected, [10, 11, 12])


class TestNewParticlesDetectionIgnoresProcessedOnes(unittest.TestCase):

    def testBoundarySecondAloneIsNotReportedAsNewWork(self):
        """Widening the watermark must not make every poll look busy."""
        clazz = _Class(1, [_Image(10, '2026-01-01 10:00:00'),
                           _Image(11, '2026-01-01 10:00:00')])
        harness = _Harness([clazz])

        harness.extract()

        self.assertFalse(
            harness._newParticlesToProcess(),
            "Every particle has been collected already: reporting new "
            "work here queues an extraction step on every single poll.",
        )

    def testGenuinelyNewParticlesAreStillReported(self):
        clazz = _Class(1, [_Image(10, '2026-01-01 10:00:00')])
        harness = _Harness([clazz])

        harness.extract()
        clazz.images.append(_Image(11, '2026-01-01 10:00:00'))

        self.assertTrue(
            harness._newParticlesToProcess(),
            "A particle sharing the boundary second is still new work.",
        )


class TestFreshOutputIsNeverReadBack(unittest.TestCase):

    def testAnUnpublishedOutputSetIsNotQueried(self):
        """A Set just created by _loadOutputSet has no storage behind it.

        Reading it back raises rather than returning an empty set, so the
        reconciliation must only touch an output already published as a
        protocol attribute.
        """
        class _UnwrittenOutput(_OutputStub):
            def getIdSet(self):
                raise AssertionError(
                    "A freshly created output Set has nothing behind it "
                    "yet: reading it back blows the whole step up."
                )

        clazz = _Class(1, [_Image(10, '2026-01-01 10:00:00')])
        harness = _Harness([clazz])
        harness._outputs = {
            OUTPUT_PARTICLES: _UnwrittenOutput(),
            OUTPUT_DISCARDED_PARTICLES: _UnwrittenOutput(),
        }

        harness.extract()

        collected = sorted(
            i.getObjId() for i in harness._outputs[OUTPUT_PARTICLES].items
        )
        self.assertEqual(collected, [10])

    def testAPublishedOutputIsReconciledAgainstOnce(self):
        """In-memory tracking can be behind what really got persisted."""
        clazz = _Class(1, [_Image(10, '2026-01-01 10:00:00')])
        harness = _Harness([clazz])
        # The output already holds this particle, but nothing in memory
        # says so - the state a resumed run starts from.
        published = harness._outputs[OUTPUT_PARTICLES]
        published.items.append(_Image(10, '2026-01-01 10:00:00'))
        setattr(harness, OUTPUT_PARTICLES, published)

        harness.extract()

        collected = sorted(i.getObjId() for i in published.items)
        self.assertEqual(
            collected,
            [10],
            "The particle is already in the output: appending it again "
            "would duplicate it in the published Set.",
        )


class TestExtractionDoesNotRescanItsOwnOutput(unittest.TestCase):

    def testPersistedIdsAreNotRebuiltFromTheOutputSetsEveryStep(self):
        clazz = _Class(1, [_Image(10, '2026-01-01 10:00:00')])
        harness = _Harness([clazz])

        harness.extract()
        baseline = harness._outputs[OUTPUT_PARTICLES].idSetCalls

        for extra in range(11, 16):
            clazz.images.append(_Image(extra, '2026-01-01 10:00:0%d'
                                       % (extra - 10)))
            harness.extract()

        self.assertEqual(
            harness._outputs[OUTPUT_PARTICLES].idSetCalls,
            baseline,
            "Rebuilding the whole persisted ID set on every extraction is "
            "O(output) per step: with millions of particles each step "
            "would re-read everything collected so far.",
        )


if __name__ == '__main__':
    unittest.main()
