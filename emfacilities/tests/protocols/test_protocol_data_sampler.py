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
from pyworkflow.object import Pointer
from emfacilities.protocols.protocol_data_sampler import ProtDataSampler


class TestDataSampler(BaseTest):
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

    def testDataSampler25(self):
        prot = self._runDataSampler("Random sampling of 0.25 proportion", batch=1000, proportion=0.25)
        self.assertSetSize(prot.outputSet, size=1175)

    def testDataSampler50(self):
        prot = self._runDataSampler("Random sampling of 0.5 proportion", batch=1000,  proportion=0.5)
        self.assertSetSize(prot.outputSet, size=2350)

    def _runDataSampler(cls, label, batch, proportion):
        protDataSampler = cls.newProtocol(ProtDataSampler,
                                          batchSize=batch,
                                          samplingProportion=proportion,
                                          delay=3)
        protDataSampler.inputImages = Pointer(cls.protImport, extended='outputParticles')
        protDataSampler.setObjLabel(label)
        cls.launchProtocol(protDataSampler)

        return protDataSampler