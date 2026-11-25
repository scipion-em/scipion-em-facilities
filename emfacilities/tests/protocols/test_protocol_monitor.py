from unittest import TestCase
from unittest.mock import patch

from emfacilities.protocols.protocol_monitor import Monitor


class TestMonitorLoop(TestCase):

    def testLoopUsesOriginalStartTimeOnContinue(self):
        class CountingMonitor(Monitor):
            def __init__(self):
                super().__init__(
                    workingDir=".",
                    samplingInterval=1,
                    monitorTime=1,
                )
                self.calls = 0

            def initLoop(self):
                pass

            def step(self):
                self.calls += 1
                return self.calls >= 2

        monitor = CountingMonitor()

        # Original monitor start: t=0. At t=90 seconds a one-minute
        # monitor has already expired. Continue must not grant it
        # another full minute.
        with patch(
                "emfacilities.protocols.protocol_monitor.time.time",
                return_value=90,
        ), patch(
                "emfacilities.protocols.protocol_monitor.time.sleep",
                return_value=None,
        ):
            monitor.loop(startTime=0)

        self.assertEqual(monitor.calls, 1)

    def testCtfMonitorUsesProtocolStartTime(self):
        from emfacilities.protocols.protocol_monitor_ctf import ProtMonitorCTF

        expectedStart = object()
        captured = {}

        class InitTime:
            def datetime(self):
                return expectedStart

        class MonitorSpy:
            def loop(self, **kwargs):
                captured.update(kwargs)

        class ProtocolSpy:
            initTime = InitTime()

            def createMonitor(self):
                return MonitorSpy()

        ProtMonitorCTF.monitorStep(ProtocolSpy())

        self.assertIs(captured.get("startTime"), expectedStart)


