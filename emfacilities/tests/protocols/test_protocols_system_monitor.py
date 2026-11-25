# ***************************************************************************
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
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

import os.path
import tempfile

import pyworkflow.tests as pwtests
import pyworkflow.utils as pwutils

import pwem.protocols as emprot
import emfacilities.protocols as monitorsProt
from emfacilities.protocols.protocol_monitor_system import MonitorSystem


class TestSystemMonitorResume(pwtests.BaseTest):

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
                        "cpu": 10.0,
                    }
                ]

        monitor = object.__new__(MonitorSystem)
        monitor._tableName = "log"
        monitor.workingDir = "/tmp"
        monitor._dataBase = "system_log.sqlite"
        monitor.timeZone = "UTC"
        monitor.timeDelta = 0
        monitor.cur = FailingCursor()

        result = monitor.getDataInflux(lastId=7)

        self.assertEqual(
            result,
            [],
            "A failed system query must not return stale cursor rows.",
        )
        self.assertFalse(
            monitor.cur.fetchallCalled,
            "fetchall() must not run after the SELECT failed.",
        )


    def testInitLoopRestoresAlertThresholdsFromExistingDatabase(self):
        with tempfile.TemporaryDirectory() as tmpDir:
            monitor = MonitorSystem(
                [],
                workingDir=tmpDir,
                samplingInterval=1,
                monitorTime=1,
                cpuAlert=50,
                memAlert=60,
                swapAlert=70,
                doGpu=False,
                doNetwork=False,
                doDiskIO=False,
                gpusToUse="0",
                nif=None,
            )
            monitor.initLoop()
            monitor.cur.execute(
                "INSERT INTO log (cpu, mem, swap) VALUES (?, ?, ?)",
                (65, 75, 72),
            )
            monitor.cur.execute(
                "INSERT INTO log (cpu, mem, swap) VALUES (?, ?, ?)",
                (80, 70, 90),
            )
            monitor.conn.close()

            resumed = MonitorSystem(
                [],
                workingDir=tmpDir,
                samplingInterval=1,
                monitorTime=1,
                cpuAlert=50,
                memAlert=60,
                swapAlert=70,
                doGpu=False,
                doNetwork=False,
                doDiskIO=False,
                gpusToUse="0",
                nif=None,
            )
            resumed.initLoop()

            try:
                self.assertEqual(resumed.cpuAlert, 80)
                self.assertEqual(resumed.memAlert, 75)
                self.assertEqual(resumed.swapAlert, 90)
            finally:
                resumed.conn.close()


class TestSystemMonitorSamplingFailureRegression(pwtests.BaseTest):
    def testNetworkSamplingFailureDoesNotCrashMonitor(self):
        from unittest.mock import patch

        class Memory:
            percent = 20.0

        with tempfile.TemporaryDirectory() as tmpDir:
            monitor = MonitorSystem(
                [],
                workingDir=tmpDir,
                samplingInterval=1,
                monitorTime=1,
                cpuAlert=101,
                memAlert=101,
                swapAlert=101,
                doGpu=False,
                doNetwork=True,
                doDiskIO=False,
                gpusToUse="0",
                nif="eth0",
            )

            try:
                with patch(
                    "emfacilities.protocols.protocol_monitor_system."
                    "psutil.cpu_percent",
                    return_value=10.0,
                ), patch(
                    "emfacilities.protocols.protocol_monitor_system."
                    "psutil.virtual_memory",
                    return_value=Memory(),
                ), patch(
                    "emfacilities.protocols.protocol_monitor_system."
                    "psutil.swap_memory",
                    return_value=Memory(),
                ), patch(
                    "emfacilities.protocols.protocol_monitor_system."
                    "psutil.net_io_counters",
                    side_effect=KeyError("eth0"),
                ):
                    monitor.initLoop()
                    finished = monitor.step()

                self.assertTrue(finished)

                monitor.cur.execute(
                    "SELECT cpu, mem, swap, eth0_send, eth0_recv FROM log"
                )
                rows = monitor.cur.fetchall()

                self.assertEqual(len(rows), 1)
                self.assertEqual(rows[0][0:3], (10.0, 20.0, 20.0))
                self.assertEqual(rows[0][3:5], (0.0, 0.0))
            finally:
                monitor.conn.close()

class TestSystemMonitorDiskSamplingFailureRegression(pwtests.BaseTest):
    def testDiskSamplingFailureDoesNotCrashMonitor(self):
        from unittest.mock import patch

        class Memory:
            percent = 20.0

        with tempfile.TemporaryDirectory() as tmpDir:
            monitor = MonitorSystem(
                [],
                workingDir=tmpDir,
                samplingInterval=1,
                monitorTime=1,
                cpuAlert=101,
                memAlert=101,
                swapAlert=101,
                doGpu=False,
                doNetwork=False,
                doDiskIO=True,
                gpusToUse="0",
                nif=None,
            )

            try:
                with patch(
                    "emfacilities.protocols.protocol_monitor_system."
                    "psutil.cpu_percent",
                    return_value=10.0,
                ), patch(
                    "emfacilities.protocols.protocol_monitor_system."
                    "psutil.virtual_memory",
                    return_value=Memory(),
                ), patch(
                    "emfacilities.protocols.protocol_monitor_system."
                    "psutil.swap_memory",
                    return_value=Memory(),
                ), patch(
                    "emfacilities.protocols.protocol_monitor_system."
                    "psutil.disk_io_counters",
                    side_effect=RuntimeError("disk counters unavailable"),
                ):
                    monitor.initLoop()
                    finished = monitor.step()

                self.assertTrue(finished)

                monitor.cur.execute(
                    "SELECT cpu, mem, swap, disk_read, disk_write FROM log"
                )
                rows = monitor.cur.fetchall()

                self.assertEqual(len(rows), 1)
                self.assertEqual(rows[0][0:3], (10.0, 20.0, 20.0))
                self.assertEqual(rows[0][3:5], (0.0, 0.0))
            finally:
                monitor.conn.close()


class TestSystemMonitorLifecycleRegression(pwtests.BaseTest):
    def testScheduledProducerKeepsSystemMonitorAlive(self):
        from unittest.mock import patch
        from pyworkflow.protocol.constants import STATUS_SCHEDULED

        class ScheduledProtocol:
            def getStatus(self):
                return STATUS_SCHEDULED

        class Memory:
            percent = 20.0

        protocol = ScheduledProtocol()

        with tempfile.TemporaryDirectory() as tmpDir:
            monitor = MonitorSystem(
                [protocol],
                workingDir=tmpDir,
                samplingInterval=1,
                monitorTime=1,
                cpuAlert=101,
                memAlert=101,
                swapAlert=101,
                doGpu=False,
                doNetwork=False,
                doDiskIO=False,
                gpusToUse="0",
                nif=None,
            )

            try:
                with patch(
                    "emfacilities.protocols.protocol_monitor_system."
                    "getUpdatedProtocol",
                    return_value=protocol,
                ), patch(
                    "emfacilities.protocols.protocol_monitor_system."
                    "psutil.cpu_percent",
                    return_value=10.0,
                ), patch(
                    "emfacilities.protocols.protocol_monitor_system."
                    "psutil.virtual_memory",
                    return_value=Memory(),
                ), patch(
                    "emfacilities.protocols.protocol_monitor_system."
                    "psutil.swap_memory",
                    return_value=Memory(),
                ):
                    monitor.initLoop()
                    finished = monitor.step()

                self.assertFalse(
                    finished,
                    "A scheduled producer is active and must keep the "
                    "system monitor alive.",
                )
            finally:
                monitor.conn.close()
