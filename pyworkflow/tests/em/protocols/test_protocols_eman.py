# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import unittest, sys
from pyworkflow.tests import *
from pyworkflow.em import *
from pyworkflow.em.packages.eman2 import *


class TestEmanBase(BaseTest):
    @classmethod
    def setData(cls, projectData='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(projectData)
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.crdsDir = cls.dataset.getFile('boxingDir')
    
    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage, scannedPixelSize, magnification, sphericalAberration):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or the ScannedPixelSize + microscope magnification
        if not samplingRate is None:
            cls.protImport = ProtImportMicrographs(samplingRateMode=0, pattern=pattern, samplingRate=samplingRate, magnification=magnification, 
                                                   voltage=voltage, sphericalAberration=sphericalAberration)
        else:
            cls.protImport = ProtImportMicrographs(samplingRateMode=1, pattern=pattern, scannedPixelSize=scannedPixelSize, 
                                                   voltage=voltage, magnification=magnification, sphericalAberration=sphericalAberration)
            
        cls.proj.launchProtocol(cls.protImport, wait=True)
        # check that input micrographs have been imported (a better way to do this?)
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. outputMicrographs is None.' % pattern)
        return cls.protImport
    
    @classmethod
    def runImportMicrographBPV(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportMicrograph(pattern, samplingRate=1.237, voltage=300, sphericalAberration=2, scannedPixelSize=None, magnification=56000)


class TestEmanBoxing(TestEmanBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanBase.setData()
        cls.protImport = cls.runImportMicrographBPV(cls.micsFn)
    
    def testCreateOutput(self):
        print "Running Eman fake particle picking..."
        protPP = EmanProtBoxing(importFolder=self.crdsDir, runMode=1)
        protPP.inputMicrographs.set(self.protImport.outputMicrographs)
        protPP.boxSize.set(550)
        self.proj.launchProtocol(protPP, wait=True)
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestEmanBoxing)
    unittest.TextTestRunner(verbosity=2).run(suite)