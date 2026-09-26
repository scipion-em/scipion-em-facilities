import os
import tempfile
from unittest import TestCase

from pyworkflow.protocol.constants import STATUS_RUNNING

from emfacilities.protocols.protocol_monitor_movie_gain import MonitorMovieGain


class DummyMovieGainProtocol:
    def __init__(self, root):
        self.root = root

    def _getPath(self, name):
        return os.path.join(self.root, name)

    def getStatus(self):
        return STATUS_RUNNING


class TestMonitorMovieGain(TestCase):

    def _newMonitor(self, root):
        return MonitorMovieGain(
            DummyMovieGainProtocol(root),
            workingDir=root,
            samplingInterval=1,
            monitorTime=1,
            stddevValue=0.04,
            ratio1Value=99,
            ratio2Value=99,
        )

    def testStepDoesNotRepeatWarningForSameSummaryLine(self):
        with tempfile.TemporaryDirectory() as tmpDir:
            summary = os.path.join(tmpDir, "summaryForMonitor.txt")
            with open(summary, "w") as handle:
                handle.write(
                    "movie_000001_residual: 0.050000 1.0 1.0 1.0\n"
                )

            monitor = self._newMonitor(tmpDir)
            monitor.initLoop()
            monitor.step()
            monitor.step()

            warnings = os.path.join(tmpDir, "warningsMonitor.txt")
            with open(warnings, "r") as handle:
                lines = handle.readlines()

            self.assertEqual(len(lines), 1)

    def testContinueRestoresProcessedSummaryLines(self):
        with tempfile.TemporaryDirectory() as tmpDir:
            summary = os.path.join(tmpDir, "summaryForMonitor.txt")
            with open(summary, "w") as handle:
                handle.write(
                    "movie_000001_residual: 0.050000 1.0 1.0 1.0\n"
                )

            first = self._newMonitor(tmpDir)
            first.initLoop()
            first.step()

            resumed = self._newMonitor(tmpDir)
            resumed.initLoop()
            resumed.step()

            warnings = os.path.join(tmpDir, "warningsMonitor.txt")
            with open(warnings, "r") as handle:
                lines = handle.readlines()

            self.assertEqual(len(lines), 1)

class TestMonitorMovieGainPartialSummaryRegression(TestCase):
    def testPartialLastSummaryLineIsRetriedWhenCompleted(self):
        with tempfile.TemporaryDirectory() as tmpDir:
            summary = os.path.join(tmpDir, "summaryForMonitor.txt")
            with open(summary, "w") as handle:
                handle.write(
                    "movie_000001_residual: 0.050000 1.0 1.0 1.0\n"
                )
                handle.write(
                    "movie_000002_residual: 0.060000 1.0"
                )

            monitor = MonitorMovieGain(
                DummyMovieGainProtocol(tmpDir),
                workingDir=tmpDir,
                samplingInterval=1,
                monitorTime=1,
                stddevValue=0.04,
                ratio1Value=99,
                ratio2Value=99,
            )
            monitor.initLoop()

            monitor.step()

            self.assertEqual(
                monitor._lastSummaryLine,
                1,
                "A partial trailing record must remain pending.",
            )

            with open(summary, "a") as handle:
                handle.write(" 1.0 1.0\n")

            monitor.step()

            self.assertEqual(monitor._lastSummaryLine, 2)

            warnings = os.path.join(tmpDir, "warningsMonitor.txt")
            with open(warnings, "r") as handle:
                lines = handle.readlines()

            self.assertEqual(len(lines), 2)

    def testPartialLastSummaryLineKeepsMonitorAliveAfterProducerFinishes(self):
        class FinishedMovieGainProtocol(DummyMovieGainProtocol):
            def getStatus(self):
                return 999999

        with tempfile.TemporaryDirectory() as tmpDir:
            summary = os.path.join(tmpDir, "summaryForMonitor.txt")
            with open(summary, "w") as handle:
                handle.write(
                    "movie_000001_residual: 0.050000 1.0"
                )

            monitor = MonitorMovieGain(
                FinishedMovieGainProtocol(tmpDir),
                workingDir=tmpDir,
                samplingInterval=1,
                monitorTime=1,
                stddevValue=0.04,
                ratio1Value=99,
                ratio2Value=99,
            )
            monitor.initLoop()

            finished = monitor.step()

            self.assertFalse(
                finished,
                "An unterminated trailing record must keep the monitor alive "
                "until that record can be retried.",
            )
            self.assertEqual(monitor._lastSummaryLine, 0)
