# ***************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# ***************************************************************************/

import os
import tempfile


import pyworkflow.tests as pwtests
import pyworkflow.protocol as pwprot

import pwem.protocols as emprot

import emfacilities.protocols as monitorsProt

# Load the number of movies for the simulation, by default equal 3,
# but can be modified in the environement
MICS = os.environ.get('SCIPION_TEST_MICS', 3)


class TestCtfStream(pwtests.BaseTest):

    def testInfluxReadFailureDoesNotReuseCursorResults(self):
        class FailingCursor:
            def __init__(self):
                self.fetchallCalled = False

            def execute(self, command):
                raise RuntimeError("simulated sqlite read failure")

            def fetchall(self):
                self.fetchallCalled = True
                return [
                    {
                        "id": 99,
                        "timestamp": "2026-09-20 10:00:00",
                    }
                ]

        monitor = object.__new__(monitorsProt.MonitorCTF)
        monitor._tableName = "log"
        monitor.workingDir = "/tmp"
        monitor._dataBase = "ctf_log.sqlite"
        monitor.timeZone = "UTC"
        monitor.timeDelta = 0
        monitor.cur = FailingCursor()

        result = monitor.getDataInflux(lastId=7)

        self.assertEqual(
            result,
            [],
            "A failed CTF query must not return stale cursor rows.",
        )
        self.assertFalse(
            monitor.cur.fetchallCalled,
            "fetchall() must not run after the SELECT failed.",
        )



    def testFailedCtfInsertIsRetriedInsteadOfMarkedAsRead(self):
        from unittest.mock import patch

        class DummyMicrograph:
            def getFileName(self):
                return "/tmp/mic.mrc"

        class DummyCtf:
            def getDefocusU(self):
                return 2000.0

            def getDefocusV(self):
                return 1500.0

            def getDefocusAngle(self):
                return 0.0

            def getResolution(self):
                return 3.0

            def getFitQuality(self):
                return 1.0

            def hasPhaseShift(self):
                return False

            def getPsdFile(self):
                return "/tmp/psd.psd"

            def getMicrograph(self):
                return DummyMicrograph()

            def getObjCreation(self):
                return "2026-09-20 10:00:00"

        class DummyCtfSet:
            def getIdSet(self):
                return {7}

            def __getitem__(self, objId):
                return DummyCtf()

        class DummyProtocol:
            outputCTF = DummyCtfSet()

            def getStatus(self):
                return 0

        class FailingCursor:
            def execute(self, sql):
                raise RuntimeError("simulated sqlite failure")

        with tempfile.TemporaryDirectory() as tmpDir:
            monitor = monitorsProt.MonitorCTF(
                DummyProtocol(),
                workingDir=tmpDir,
                samplingInterval=1,
                monitorTime=1,
                minDefocus=1000,
                maxDefocus=40000,
                astigmatism=2000,
            )
            monitor.initLoop()
            monitor.cur = FailingCursor()

            try:
                with patch(
                    "emfacilities.protocols.protocol_monitor_ctf.getUpdatedProtocol",
                    return_value=DummyProtocol(),
                ):
                    monitor.step()

                self.assertNotIn(
                    7,
                    monitor.readCTFs,
                    "A CTF whose log insert failed must remain pending for retry.",
                )
            finally:
                monitor.conn.close()



    def testInitLoopRestoresReadCtfIdsFromExistingDatabase(self):
        class DummyProtocol:
            pass

        with tempfile.TemporaryDirectory() as tmpDir:
            monitor = monitorsProt.MonitorCTF(
                DummyProtocol(),
                workingDir=tmpDir,
                samplingInterval=1,
                monitorTime=1,
                minDefocus=1000,
                maxDefocus=40000,
                astigmatism=2000,
            )
            monitor.initLoop()
            monitor.cur.execute(
                "INSERT INTO log (ctfID) VALUES (?)",
                (7,),
            )
            monitor.conn.close()

            resumed = monitorsProt.MonitorCTF(
                DummyProtocol(),
                workingDir=tmpDir,
                samplingInterval=1,
                monitorTime=1,
                minDefocus=1000,
                maxDefocus=40000,
                astigmatism=2000,
            )
            resumed.initLoop()

            try:
                self.assertEqual(resumed.readCTFs, {7})
            finally:
                resumed.conn.close()



    def testInitLoopRestoresDefocusAlertThresholdsFromExistingDatabase(self):
        class DummyProtocol:
            pass

        with tempfile.TemporaryDirectory() as tmpDir:
            monitor = monitorsProt.MonitorCTF(
                DummyProtocol(),
                workingDir=tmpDir,
                samplingInterval=1,
                monitorTime=1,
                minDefocus=1000,
                maxDefocus=40000,
                astigmatism=2000,
            )
            monitor.initLoop()
            monitor.cur.execute(
                "INSERT INTO log (ctfID, defocusU, defocusV) VALUES (?, ?, ?)",
                (7, 45000, 900),
            )
            monitor.cur.execute(
                "INSERT INTO log (ctfID, defocusU, defocusV) VALUES (?, ?, ?)",
                (8, 47000, 850),
            )
            monitor.conn.close()

            resumed = monitorsProt.MonitorCTF(
                DummyProtocol(),
                workingDir=tmpDir,
                samplingInterval=1,
                monitorTime=1,
                minDefocus=1000,
                maxDefocus=40000,
                astigmatism=2000,
            )
            resumed.initLoop()

            try:
                self.assertEqual(resumed.maxDefocus, 47000)
                self.assertEqual(resumed.minDefocus, 850)
            finally:
                resumed.conn.close()


    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)

    def _updateProtocol(self, prot):
        prot2 = pwprot.getProtocolFromDb(prot.getProject().path,
                                         prot.getDbPath(),
                                         prot.getObjId())
        # Close DB connections
        prot2.getProject().closeMapper()
        prot2.closeMappers()
        return prot2

    def test_pattern(self):
        """ Import several Particles from a given pattern.
        """
        kwargs = {'xDim': 4096,
                  'yDim': 4096,
                  'nDim': MICS,
                  'samplingRate': 1.25,
                  'creationInterval': 5,
                  'delay': 0,
                  'setof': emprot.SET_OF_RANDOM_MICROGRAPHS  # SetOfMicrographs
                  }

        # put some stress on the system
        protStream = self.newProtocol(emprot.ProtCreateStreamData, **kwargs)
        protStream.setObjLabel('create Stream Mic')
        self.proj.launchProtocol(protStream, wait=False)

        self._waitOutput(protStream, 'outputMicrographs')

        # then introduce monitor, checking all the time ctf and saving to
        # database
        kwargs = {
            'useCtffind4': True,
            'ctfDownFactor': 2,
            'numberOfThreads': 4
        }
        from pwem import Domain

        ProtCTFFind = Domain.importFromPlugin('cistem.protocols',
                                              'CistemProtCTFFind', doRaise=True)
        protCTF = self.newProtocol(ProtCTFFind, **kwargs)
        protCTF.inputMicrographs.set(protStream.outputMicrographs)
        self.proj.launchProtocol(protCTF, wait=False)

        self._waitOutput(protCTF, 'outputCTF')

        kwargs = {'samplingInterval': 10,
                  'interval': 300,
                  'maxDefocus': 40000,
                  'minDefocus': 1000,
                  'astigmatism': 0.2,
                  'monitorTime': 5
                  }

        protMonitor = self.newProtocol(monitorsProt.ProtMonitorCTF, **kwargs)
        protMonitor.inputProtocol.set(protCTF)
        self.launchProtocol(protMonitor)

        baseFn = protMonitor._getPath(monitorsProt.CTF_LOG_SQLITE)
        self.assertTrue(os.path.isfile(baseFn))
