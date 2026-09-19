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
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from datetime import datetime
from unittest.mock import patch
from pwem.protocols.protocol_import import ProtImportParticles
import pwem.protocols as emprot
from pyworkflow.object import Pointer
from pyworkflow.object import Set
from emfacilities.protocols.protocol_good_classes_extractor import ProtGoodClassesExtractor


class TestGoodClassesExtractor(BaseTest):
    """ Test good classes extractor protocol """

    @classmethod
    def setData(cls):
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.setData()
        # Run needed protocols
        cls.runImportParticles()
        cls.runClassSelector()

    @classmethod
    def runImportParticles(cls):
        """
        Import an EMX file with Particles and defocus
        """
        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         objLabel='from relion (classify 2d)',
                                         importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                         starFile=cls.dsRelion.getFile('import/classify2d/extra/relion_it015_data.star'),
                                         magnification=10000,
                                         samplingRate=7.08,
                                         haveDataBeenPhaseFlipped=True
                                         )

        cls.launchProtocol(cls.protImport)

    @classmethod
    def runClassSelector(cls):
        """
        Add tests for classes selector representatives
        """
        cls.classSelector = cls.newProtocol(emprot.ProtClassesSelector,
                                            objLabel='representatives from 10 mayor classes',
                                            firstNElements=10,
                                            extractRepresentative=True)

        cls.classSelector.inputClasses = Pointer(cls.protImport, extended='outputClasses')
        cls.launchProtocol(cls.classSelector)




    def testContinueDetectsClosedStreamWithoutNewParticles(self):
        class ClosedInputSet:
            def __init__(self):
                self.closed = False

            def getFileName(self):
                return "classes.sqlite"

            def getStreamState(self):
                return Set.STREAM_CLOSED

            def close(self):
                self.closed = True

        class InputPointer:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        inputSet = ClosedInputSet()

        prot = self.newProtocol(ProtGoodClassesExtractor)
        prot.inputClasses = InputPointer(inputSet)
        prot.dictsTimes = {"1": "2026-09-19 08:00:00"}
        prot.isStreamClosed = Set.STREAM_OPEN
        prot.lastCheck = datetime.now()

        prot._loadInputClassesSet = lambda: inputSet

        with patch(
                "emfacilities.protocols.protocol_good_classes_extractor.os.path.getmtime",
                return_value=0,
        ):
            hasNewParticles = prot._newParticlesToProcess()

        self.assertFalse(hasNewParticles)
        self.assertEqual(prot.isStreamClosed, Set.STREAM_CLOSED)
        self.assertTrue(inputSet.closed)


    def testContinueRestoresLastProcessedClassTimes(self):
        expectedTimes = {
            "1": "2026-09-19 08:00:00",
            "5": "2026-09-19 08:01:00",
        }

        prot = self.newProtocol(ProtGoodClassesExtractor)
        prot.isContinued = lambda: True
        prot._getLastDone = lambda: expectedTimes.copy()

        prot.initialStep()

        self.assertEqual(prot.dictsTimes, expectedTimes)


    def testGoodClassesSelectorAvgs(self):
        prot = self._runGoodClassesSelectorAverages("Select good particles from averages")
        self.assertSetSize(prot.outputParticles, size=4165)
        self.assertSetSize(prot.outputParticlesDiscarded, size=535)


    def testGoodClassesSelectorIds(self):
        prot = self._runGoodClassesSelectorIds("Select good particles from list ids")
        self.assertSetSize(prot.outputParticles, size=4165)
        self.assertSetSize(prot.outputParticlesDiscarded, size=535)

    def _runGoodClassesSelectorAverages(cls, label):
        protGoodClassSelectorAvg = cls.newProtocol(ProtGoodClassesExtractor)

        protGoodClassSelectorAvg.inputClasses = Pointer(cls.protImport, extended='outputClasses')
        protGoodClassSelectorAvg.inputGoodClasses = Pointer(cls.classSelector, extended='output')
        protGoodClassSelectorAvg.setObjLabel(label)
        cls.launchProtocol(protGoodClassSelectorAvg)

        return protGoodClassSelectorAvg

    def _runGoodClassesSelectorIds(cls, label):
        protGoodClassSelectorIds = cls.newProtocol(ProtGoodClassesExtractor,
                                                   mode=ProtGoodClassesExtractor.LIST_IDS,
                                                   inputGoodListIds="1,5,6,16,17,18,20,24,25,31")

        protGoodClassSelectorIds.inputClasses = Pointer(cls.protImport, extended='outputClasses')
        protGoodClassSelectorIds.setObjLabel(label)
        cls.launchProtocol(protGoodClassSelectorIds)

        return protGoodClassSelectorIds
    def testContinueRestoresParticleCountersFromOutputs(self):
        class OutputSet:
            def __init__(self, ids):
                self._ids = set(ids)

            def getIdSet(self):
                return set(self._ids)

        prot = self.newProtocol(ProtGoodClassesExtractor)
        prot.isContinued = lambda: True
        prot._getLastDone = lambda: {
            "1": "2026-09-19 08:00:00",
        }
        prot.outputParticles = OutputSet({1, 2, 3})
        prot.outputParticlesDiscarded = OutputSet({4, 5})

        prot.initialStep()

        self.assertEqual(set(prot.goodParticles), {1, 2, 3})
        self.assertEqual(set(prot.badParticles), {4, 5})
        self.assertEqual(
            prot.particlesDistribution,
            {"good": [3], "bad": [2]},
        )

    def testContinueClassWithoutNewParticlesKeepsCheckpoint(self):
        class EmptyOutput:
            def __len__(self):
                return 0

            def append(self, item):
                raise AssertionError("No particle should be appended")

        class EmptyClass:
            def getObjId(self):
                return 1

            def iterItems(self, **kwargs):
                return iter(())

        class InputClasses:
            def iterItems(self, **kwargs):
                return iter((EmptyClass(),))

        prot = self.newProtocol(ProtGoodClassesExtractor)
        previousTime = "2026-09-19 08:00:00"
        prot.dictsTimes = {"1": previousTime}
        prot.goodClassesIDs = [1]
        prot.goodParticles = []
        prot.badParticles = []
        prot.particlesDistribution = {"good": [], "bad": []}
        prot.isStreamClosed = Set.STREAM_OPEN

        prot._loadOutputSet = lambda *args: EmptyOutput()
        prot._updateOutputSet = lambda *args, **kwargs: None
        prot._writeLastDone = lambda value: None
        prot._createPlots = lambda: None

        prot.extractElements(InputClasses())

        self.assertEqual(prot.dictsTimes["1"], previousTime)
        self.assertEqual(prot.goodParticles, [])
        self.assertEqual(prot.badParticles, [])




