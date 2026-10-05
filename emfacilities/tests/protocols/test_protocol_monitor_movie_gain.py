import os
import tempfile
from unittest import TestCase
from unittest.mock import patch

from pyworkflow.protocol.constants import STATUS_RUNNING

from emfacilities.protocols.protocol_monitor_movie_gain import MonitorMovieGain


class MovieGainMonitorTestCase(TestCase):
    def setUp(self):
        self._updatedProtocolPatcher = patch(
            "emfacilities.protocols.protocol_monitor_movie_gain."
            "getUpdatedProtocol",
            side_effect=lambda protocol: protocol,
        )
        self._updatedProtocolPatcher.start()
        self.addCleanup(self._updatedProtocolPatcher.stop)

class DummyMovieGainProtocol:
    def __init__(self, root):
        self.root = root

    def _getPath(self, name):
        return os.path.join(self.root, name)

    def getStatus(self):
        return STATUS_RUNNING


class TestMonitorMovieGain(MovieGainMonitorTestCase):

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

    def testStepReadsOnlyNewSummaryTailAfterFirstPoll(self):
        from unittest.mock import patch

        with tempfile.TemporaryDirectory() as tmpDir:
            summary = os.path.join(tmpDir, "summaryForMonitor.txt")
            with open(summary, "w") as handle:
                handle.write(
                    "movie_000001_residual: 0.010000 1.0 1.0 1.0\n"
                )
                handle.write(
                    "movie_000002_residual: 0.010000 1.0 1.0 1.0\n"
                )

            monitor = self._newMonitor(tmpDir)
            monitor.initLoop()
            monitor.step()

            with open(summary, "a") as handle:
                handle.write(
                    "movie_000003_residual: 0.010000 1.0 1.0 1.0\n"
                )

            realOpen = open
            seekPositions = []
            readlineCalls = []

            class _IncrementalReader:
                def __init__(self, handle):
                    self._handle = handle

                def __enter__(self):
                    self._handle.__enter__()
                    return self

                def __exit__(self, *args):
                    return self._handle.__exit__(*args)

                def seek(self, offset, whence=0):
                    seekPositions.append((offset, whence))
                    return self._handle.seek(offset, whence)

                def tell(self):
                    return self._handle.tell()

                def readline(self, *args, **kwargs):
                    readlineCalls.append(True)
                    return self._handle.readline(*args, **kwargs)

                def readlines(self, *args, **kwargs):
                    raise AssertionError(
                        "Movie gain polling must not reread the full "
                        "summary file on every step."
                    )

                def read(self, *args, **kwargs):
                    raise AssertionError(
                        "Movie gain polling must read incrementally by line."
                    )

                def __iter__(self):
                    raise AssertionError(
                        "Movie gain polling must seek to the saved offset "
                        "instead of iterating from the beginning."
                    )

            def incrementalOpen(path, mode="r", *args, **kwargs):
                handle = realOpen(path, mode, *args, **kwargs)
                if (
                    os.path.abspath(path) == os.path.abspath(summary)
                    and "r" in mode
                ):
                    return _IncrementalReader(handle)
                return handle

            with patch("builtins.open", side_effect=incrementalOpen):
                monitor.step()

            self.assertTrue(
                any(offset > 0 for offset, _ in seekPositions),
                "The second poll must seek to the previously consumed "
                "summary offset.",
            )
            self.assertLessEqual(
                len(readlineCalls),
                2,
                "Only the newly appended line and EOF should be read.",
            )



    def testScheduledProducerWithSummaryKeepsMonitorAlive(self):
        from unittest.mock import patch
        from pyworkflow.protocol.constants import STATUS_SCHEDULED

        class ScheduledMovieGainProtocol(DummyMovieGainProtocol):
            def getStatus(self):
                return STATUS_SCHEDULED

        with tempfile.TemporaryDirectory() as tmpDir:
            summary = os.path.join(tmpDir, "summaryForMonitor.txt")
            with open(summary, "w") as handle:
                handle.write(
                    "movie_000001_residual: 0.010000 1.0 1.0 1.0\n"
                )

            monitor = MonitorMovieGain(
                ScheduledMovieGainProtocol(tmpDir),
                workingDir=tmpDir,
                samplingInterval=1,
                monitorTime=1,
                stddevValue=0.04,
                ratio1Value=99,
                ratio2Value=99,
            )
            monitor.initLoop()

            with patch(
                "emfacilities.protocols.protocol_monitor_movie_gain."
                "getUpdatedProtocol",
                return_value=ScheduledMovieGainProtocol(tmpDir),
            ):
                result = monitor.step()

            self.assertFalse(
                result,
                "A scheduled movie-gain producer is still active and must "
                "not stop its monitor before it starts running.",
            )

    def testFinishedProducerWithoutSummaryStopsMonitor(self):
        from unittest.mock import patch
        class FinishedMovieGainProtocol(DummyMovieGainProtocol):
            def getStatus(self):
                return "finished"

        with tempfile.TemporaryDirectory() as tmpDir:
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

            with patch(
                "emfacilities.protocols.protocol_monitor_movie_gain."
                "getUpdatedProtocol",
                return_value=FinishedMovieGainProtocol(tmpDir),
            ):
                result = monitor.step()

            self.assertTrue(
                result,
                "A finished movie-gain producer with no summary file must "
                "not keep the monitor alive until monitorTime expires.",
            )

    def testStepRefreshesProducerStatusBeforeStopping(self):
        from unittest.mock import patch

        class FinishedMovieGainProtocol(DummyMovieGainProtocol):
            def getStatus(self):
                return "finished"

        with tempfile.TemporaryDirectory() as tmpDir:
            summary = os.path.join(tmpDir, "summaryForMonitor.txt")
            with open(summary, "w") as handle:
                handle.write(
                    "movie_000001_residual: 0.010000 1.0 1.0 1.0\n"
                )

            originalProtocol = DummyMovieGainProtocol(tmpDir)
            refreshedProtocol = FinishedMovieGainProtocol(tmpDir)

            monitor = MonitorMovieGain(
                originalProtocol,
                workingDir=tmpDir,
                samplingInterval=1,
                monitorTime=1,
                stddevValue=0.04,
                ratio1Value=99,
                ratio2Value=99,
            )
            monitor.initLoop()

            with patch(
                "emfacilities.protocols.protocol_monitor_movie_gain."
                "getUpdatedProtocol",
                return_value=refreshedProtocol,
                create=True,
            ):
                finished = monitor.step()

            self.assertTrue(
                finished,
                "Movie-gain polling must refresh the producer protocol "
                "before deciding whether monitoring is complete.",
            )

class TestMonitorMovieGainPartialSummaryRegression(MovieGainMonitorTestCase):
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


class TestMonitorMovieGainLogicalResumeRegression(TestMonitorMovieGain):
    def testResumeDoesNotDependOnExternalLastLineCheckpoint(self):
        with tempfile.TemporaryDirectory() as tmpDir:
            summary = os.path.join(tmpDir, "summaryForMonitor.txt")
            checkpoint = os.path.join(tmpDir, "movie_gain_monitor.last_line")
            with open(summary, "w") as handle:
                handle.write("movie_000001_residual: 0.050000 1.0 1.0 1.0\n")
                handle.write("movie_000002_residual: 0.060000 1.0 1.0 1.0\n")
            with open(checkpoint, "w") as handle:
                handle.write("999")

            monitor = self._newMonitor(tmpDir)
            monitor.initLoop()
            monitor.step()

            warnings = os.path.join(tmpDir, "warningsMonitor.txt")
            with open(warnings, "r") as handle:
                lines = handle.readlines()

            self.assertEqual(len(lines), 2)
            self.assertEqual(monitor._lastSummaryLine, 2)

    def testMonitorDoesNotCreateLastLineCheckpoint(self):
        with tempfile.TemporaryDirectory() as tmpDir:
            summary = os.path.join(tmpDir, "summaryForMonitor.txt")
            with open(summary, "w") as handle:
                handle.write("movie_000001_residual: 0.050000 1.0 1.0 1.0\n")

            monitor = self._newMonitor(tmpDir)
            monitor.initLoop()
            monitor.step()

            self.assertFalse(os.path.exists(os.path.join(tmpDir, "movie_gain_monitor.last_line")))
