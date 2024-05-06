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
from pwem.protocols.protocol_good_class_extractor import ProtGoodClassesExtractor


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