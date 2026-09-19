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
from datetime import datetime, timedelta
from pyworkflow.tests import BaseTest, DataSet
from pwem.protocols.protocol_import import ProtImportMicrographs
from pyworkflow.object import Pointer
import pyworkflow.tests as tests
from emfacilities.protocols.protocol_data_counter import ProtDataCounter


class TestDataCounter(BaseTest):
    """ Test data counter protocol """

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





    def testTimeoutDaysUse24Hours(self):
        prot = self.newProtocol(ProtDataCounter)

        self.assertEqual(prot.getTimeOutInSeconds("1d"), 86400)
        self.assertEqual(
            prot.getTimeOutInSeconds("1d 2h 20m 15s"),
            86400 + 2 * 3600 + 20 * 60 + 15,
        )


    def testTimerUsesPreservedProtocolStartOnContinue(self):
        class SummaryVar:
            def __init__(self):
                self.value = None

            def set(self, value):
                self.value = value

        prot = self.newProtocol(
            ProtDataCounter,
            outputSize=100,
            boolTimer=True,
            timeout="10s",
        )
        prot.finished = False
        prot.timerOut = False
        prot.timeoutSecs = 10
        prot.lastTimeCheckTimer = datetime.now()
        prot.summaryVar = SummaryVar()
        prot.initTime.set(datetime.now() - timedelta(seconds=7))

        prot.timerStep()

        self.assertFalse(prot.timerOut)
        self.assertLessEqual(
            prot.timeoutSecs,
            3,
            "Continue must preserve the elapsed timer budget from the original run.",
        )


    def testTimerExpiresWithoutNewInput(self):
        class SummaryVar:
            def __init__(self):
                self.value = None

            def set(self, value):
                self.value = value

        prot = self.newProtocol(
            ProtDataCounter,
            outputSize=100,
            boolTimer=True,
            timeout="10s",
        )
        prot.finished = False
        prot.timerOut = False
        prot.timeoutSecs = 10
        prot.lastTimeCheckTimer = datetime.now() - timedelta(seconds=11)
        prot.summaryVar = SummaryVar()

        # Simulate an idle streaming round: no new input and no new output.
        prot._checkNewInput = lambda: None
        prot._checkNewOutput = lambda: None

        prot._stepsCheck()

        self.assertTrue(
            prot.timerOut,
            "The timer must expire even when no new input batch arrives.",
        )


    def testDataCounter2000(self):
        prot = self._runDataCounter("Counter images till 1", outputSize=1)
        self.assertSetSize(prot.outputSet, size=1)


    def testDataCounter4000(self):
        prot = self._runDataCounter("Counter images till 2", outputSize=2)
        self.assertSetSize(prot.outputSet, size=2)

    def _runDataCounter(cls, label, outputSize):
        protDataSampler = cls.newProtocol(ProtDataCounter,
                                          outputSize=outputSize,
                                          delay=3)
        protDataSampler.inputImages = Pointer(cls.protImport, extended='outputMicrographs')
        protDataSampler.setObjLabel(label)
        cls.launchProtocol(protDataSampler)

        return protDataSampler