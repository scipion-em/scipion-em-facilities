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
import unittest
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols.protocol_import import ProtImportParticles
import pwem.protocols as emprot
from pyworkflow.object import Pointer
from emfacilities.protocols.protocol_volume_extractor import ProtVolumeExtractor



class TestVolumeExtractorInstanceState(unittest.TestCase):

    def testOutputsToDefineAreNotSharedAcrossProtocolInstances(self):
        class OutputParticles:
            def write(self):
                pass

        class OutputVolume:
            pass

        def newProtocolRecorder():
            protocol = object.__new__(ProtVolumeExtractor)
            recordedOutputs = {}
            object.__setattr__(
                protocol,
                "_defineOutputs",
                lambda **kwargs: recordedOutputs.update(kwargs),
            )
            object.__setattr__(
                protocol,
                "_store",
                lambda output: None,
            )
            return protocol, recordedOutputs

        ProtVolumeExtractor.outputsToDefine.clear()

        first, firstOutputs = newProtocolRecorder()
        first.createOutput(OutputParticles(), None)

        second, secondOutputs = newProtocolRecorder()
        second.createOutput(None, OutputVolume())

        self.assertEqual(
            set(firstOutputs),
            {"outputParticles"},
        )
        self.assertEqual(
            set(secondOutputs),
            {"bestVolume"},
            "Each protocol instance must define only its own outputs.",
        )


class TestVolumeExtractorValidateRegression(unittest.TestCase):
    # Regression test: Set.getItem raises rather than returning None for a
    # row it cannot find. A user-entered volume reference id that does not
    # exist in the input classes must be caught by _validate() with a clear
    # message, instead of crashing extractElements() with a confusing raw
    # exception once the protocol is launched.

    class _Value:
        def __init__(self, value):
            self._value = value

        def get(self):
            return self._value

    class _InputClasses:
        def __init__(self, ids):
            self._ids = set(ids)

        def __contains__(self, objId):
            return objId in self._ids

    def _newProtocol(self, selectBig, selectID, volumeID, inputClasses):
        protocol = object.__new__(ProtVolumeExtractor)
        object.__setattr__(protocol, "selectBig", self._Value(selectBig))
        object.__setattr__(protocol, "selectID", self._Value(selectID))
        object.__setattr__(protocol, "volumeID", self._Value(volumeID))
        object.__setattr__(protocol, "inputClasses", self._Value(inputClasses))
        return protocol

    def testValidateRejectsReferenceIdNotInInputClasses(self):
        protocol = self._newProtocol(
            selectBig=False,
            selectID=True,
            volumeID=99,
            inputClasses=self._InputClasses(ids=[1, 2, 3]),
        )

        errors = protocol._validate()

        self.assertEqual(1, len(errors))

    def testValidateAcceptsReferenceIdPresentInInputClasses(self):
        protocol = self._newProtocol(
            selectBig=False,
            selectID=True,
            volumeID=2,
            inputClasses=self._InputClasses(ids=[1, 2, 3]),
        )

        errors = protocol._validate()

        self.assertEqual([], errors)

    def testValidateSkipsCheckWhenSelectingBiggestClass(self):
        protocol = self._newProtocol(
            selectBig=True,
            selectID=False,
            volumeID=99,
            inputClasses=self._InputClasses(ids=[1, 2, 3]),
        )

        errors = protocol._validate()

        self.assertEqual([], errors)


class TestVolumeExtractor(BaseTest):
    """ Test good classes extractor protocol """

    @classmethod
    def setData(cls):
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.setData()
        # Run needed protocols
        cls.importFromRelionRefine3D()

    @classmethod
    def importFromRelionRefine3D(cls):
        """ Import aligned Particles
        """
        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         objLabel='particles from relion (auto-refine 3d)',
                                         importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                         starFile=
                                         cls.dsRelion.getFile('import/classify3d/extra/relion_it015_data.star'),
                                         magnification=10000,
                                         samplingRate=7.08,
                                         haveDataBeenPhaseFlipped=True)
        cls.launchProtocol(cls.protImport)

    def testSelectBiggestClass(self):
        prot = self._runBiggestVolumeExtractor("Select the class with the most particles")
        self.assertSetSize(prot.outputParticles, size=2989)
        self.assertIsNotNone(prot.bestVolume, "The volume does not exists")

    def testSelectIdClass(self):
        prot = self._runIdVolumeExtractor("Select the class from the ref ID")
        self.assertSetSize(prot.outputParticles, size=2989)
        self.assertIsNotNone(prot.bestVolume, "The volume does not exists")

    def _runBiggestVolumeExtractor(cls, label):
        protVol = cls.newProtocol(ProtVolumeExtractor)
        protVol.inputClasses = Pointer(cls.protImport, extended='outputClasses')
        protVol.setObjLabel(label)
        cls.launchProtocol(protVol)

        return protVol

    def _runIdVolumeExtractor(cls, label):
        protVol = cls.newProtocol(ProtVolumeExtractor,
                                  selectBig=False,
                                  selectID=True,
                                  volumeID=2
                                  )
        protVol.inputClasses = Pointer(cls.protImport, extended='outputClasses')
        protVol.setObjLabel(label)
        cls.launchProtocol(protVol)

        return protVol
