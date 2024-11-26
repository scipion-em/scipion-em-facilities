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
from pyworkflow.tests import BaseTest, DataSet
from pwem.protocols.protocol_import import ProtImportMicrographs
from pyworkflow.object import Pointer
import pyworkflow.tests as tests
from emfacilities.protocols.protocol_data_sampler import ProtDataSampler


class TestDataSampler(BaseTest):
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

    def _runDataSampler(cls, label, batch, proportion):
        protDataSampler = cls.newProtocol(ProtDataSampler,
                                          batchSize=batch,
                                          samplingProportion=proportion,
                                          delay=3)
        protDataSampler.inputImages = Pointer(cls.protImport, extended='outputMicrographs')
        protDataSampler.setObjLabel(label)
        cls.launchProtocol(protDataSampler)

        return protDataSampler