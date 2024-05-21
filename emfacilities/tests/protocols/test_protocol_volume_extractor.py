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
from pwem.protocols.protocol_import import ProtImportParticles
import pwem.protocols as emprot
from pyworkflow.object import Pointer
from emfacilities.protocols.protocol_volume_extractor import ProtVolumeExtractor


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
