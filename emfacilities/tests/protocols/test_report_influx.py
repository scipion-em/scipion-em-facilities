import os
import sys
import tempfile
import types
import unittest
from datetime import datetime
from unittest.mock import patch

from emfacilities.constants import EMFACILITIES_HOME_VARNAME
from emfacilities.protocols.report_influx import CONFILE, ReportInflux
from emfacilities.protocols.transport import Connect



_REPORT_SECRETS = """[influx]
dataBase = scipion
passwordInflux = cGFzcw==
usernameInflux = dXNlcg==
hostinflux = localhost
port = 8086
ssl = false
verify_ssl = false
TimeDelta = 0
apacheImgDir = /tmp
[paramiko]
usernameParamiko = user
passwordParamiko = pass
keyfilepath = key
keyfiletype = rsa
remote_path = /tmp
hostparamiko = localhost
"""


class _ReportProject:
    def __init__(self, shortName):
        self.shortName = shortName

    def getShortName(self):
        return self.shortName

    def getCreationTime(self):
        return datetime(2026, 1, 1, 12, 0, 0)


class _ReportProtocol:
    def __init__(self, workingDir, projectName):
        self.workingDir = workingDir
        self.reportPath = os.path.join(workingDir, "index.html")
        self.reportDir = workingDir
        self.project = _ReportProject(projectName)

    def _getCtfProtocol(self):
        return None

    def _getAlignProtocol(self):
        return None

    def _getTmpPath(self, name):
        return os.path.join(self.workingDir, name)

    def getProject(self):
        return self.project

    def getInputProtocols(self):
        return []


class _ReportProvider:
    acquisition = []

    def refreshObjects(self):
        pass

    def getObjects(self):
        return []


class _SystemMonitor:
    def getData(self, lastId=-1):
        return [{
            "id": 1,
            "timestamp": datetime(2026, 1, 1, 12, 0, 1),
            "cpu": 10.0,
        }]


class _MovieGainMonitor:
    def getData(self, lastId=-1):
        return [
            {"idx": 1, "stddev": 0.1, "ratio1": 1.1, "ratio2": 1.2},
            {"idx": 2, "stddev": 0.2, "ratio1": 1.3, "ratio2": 1.4},
        ]


class _InfluxClient:
    instance = None
    failGainWriteNumber = None

    def __init__(self, **kwargs):
        type(self).instance = self
        self.switchedDatabase = None
        self.deletedMeasurements = []
        self.measurements = []
        self.gainTimes = []
        self.gainWrites = 0

    def switch_database(self, name):
        self.switchedDatabase = name

    def delete_series(self, measurement):
        self.deletedMeasurements.append(measurement)

    def create_retention_policy(self, *args, **kwargs):
        pass

    def write_points(self, points):
        for point in points:
            self.measurements.append(point["measurement"])

            if point.get("tags", {}).get("section") != "gain":
                continue

            self.gainWrites += 1
            self.gainTimes.append(point["time"])

            if self.gainWrites == self.failGainWriteNumber:
                raise RuntimeError("simulated influx failure")


class TestReportInfluxResume(unittest.TestCase):

    def _createReport(
            self,
            tmpDir,
            projectName,
            sysMonitor=None,
            movieGainMonitor=None,
            transferFilesResult=True,
            failGainWriteNumber=None,
            existingState=None,
    ):
        with open(os.path.join(tmpDir, "secrets.cfg"), "w") as handle:
            handle.write(_REPORT_SECRETS)

        if existingState is not None:
            with open(os.path.join(tmpDir, CONFILE), "w") as handle:
                handle.write(existingState)

        influxClientClass = type(
            "FakeInfluxDBClient",
            (_InfluxClient,),
            {"failGainWriteNumber": failGainWriteNumber},
        )

        fakeInfluxModule = types.ModuleType("influxdb")
        fakeInfluxModule.InfluxDBClient = influxClientClass

        with patch.dict(
                os.environ,
                {EMFACILITIES_HOME_VARNAME: tmpDir},
        ), patch.dict(
                sys.modules,
                {"influxdb": fakeInfluxModule},
        ):
            report = ReportInflux(
                _ReportProtocol(tmpDir, projectName),
                ctfMonitor=None,
                sysMonitor=sysMonitor,
                movieGainMonitor=movieGainMonitor,
            )

        report.provider = _ReportProvider()
        report.transferFiles = lambda: transferFilesResult
        return report, influxClientClass

    @staticmethod
    def _generate(report):
        with patch.dict(
                os.environ,
                {"SCIPION_VERSION": "test"},
        ):
            return report.generate(finished=False)


    def testTransferFilesReportsPendingWorkWhenRefreshWindowExpires(self):
        class FakeResult:
            def __len__(self):
                return 1

            def get_points(self):
                return iter([
                    {
                        "time": "2026-09-20T10:00:00Z",
                        "id": 7,
                        "section": "ctf",
                        "transferImage": False,
                        "shiftPlotPathLocal": "/tmp/shift.png",
                        "psdPathLocal": "/tmp/psd.mrc",
                        "psdPathLocalPng": "/tmp/psd.png",
                        "micPathLocal": "/tmp/mic.mrc",
                        "micPathLocalPng": "/tmp/mic.png",
                        "micPath": "mic.png",
                        "psdPath": "psd.png",
                    }
                ])

        class FakeClient:
            def __init__(self):
                self.queryCalls = 0

            def query(self, query):
                self.queryCalls += 1
                return FakeResult()

            def write_points(self, points):
                pass

        class FakeImageHandler:
            def convert(self, source, target):
                pass

        class FakeConnect:
            def __init__(self, **kwargs):
                pass

            def put(self, source, target):
                pass

            def close(self):
                pass

        report = object.__new__(ReportInflux)
        report.projectName = "project"
        report.client = FakeClient()
        report.hostparamiko = "localhost"
        report.usernameParamiko = "dXNlcg=="
        report.keyfilepath = "a2V5"
        report.keyfiletype = "cnNh"
        report.remote_path = "/remote"
        report.ih = FakeImageHandler()
        report.refreshSecs = 1

        with patch(
            "emfacilities.protocols.report_influx.Connect",
            FakeConnect,
        ), patch(
            "emfacilities.protocols.report_influx.time.time",
            side_effect=[100.0, 102.0],
        ):
            self.assertFalse(
                report.transferFiles(),
                "Expiring the refresh window with pending images must not "
                "report transfer completion.",
            )


    def testGenerateDoesNotFinishWhenFileTransferFails(self):
        with tempfile.TemporaryDirectory() as tmpDir:
            report, _ = self._createReport(
                tmpDir,
                projectName="transfer-project",
                sysMonitor=_SystemMonitor(),
                transferFilesResult=False,
            )

            self.assertFalse(
                self._generate(report),
                "A failed image transfer must keep the summary monitor running.",
            )



    def testGainProgressIsPersistedAfterEachSuccessfulWrite(self):
        with tempfile.TemporaryDirectory() as tmpDir:
            report, _ = self._createReport(
                tmpDir,
                projectName="gain-resume-project",
                movieGainMonitor=_MovieGainMonitor(),
                failGainWriteNumber=2,
            )

            with self.assertRaisesRegex(
                    RuntimeError,
                    "simulated influx failure",
            ):
                self._generate(report)

            report.confParser.read(report.confFileName)
            self.assertEqual(
                report.confParser.getint("gain", "lastId"),
                1,
                "A successful gain write must be checkpointed before the next item.",
            )



    def testGainPointsUseDistinctTimestampsWithinSameBatch(self):
        with tempfile.TemporaryDirectory() as tmpDir:
            report, influxClientClass = self._createReport(
                tmpDir,
                projectName="gain-project",
                movieGainMonitor=_MovieGainMonitor(),
            )
            self._generate(report)

            gainTimes = influxClientClass.instance.gainTimes
            self.assertEqual(len(gainTimes), 2)
            self.assertNotEqual(
                gainTimes[0],
                gainTimes[1],
                "Gain points in the same batch must not share an Influx timestamp.",
            )



    def testGenerateKeepsSlugifiedMeasurementName(self):
        with tempfile.TemporaryDirectory() as tmpDir:
            report, influxClientClass = self._createReport(
                tmpDir,
                projectName="Project #1",
            )
            self._generate(report)

            self.assertTrue(influxClientClass.instance.measurements)
            self.assertEqual(
                set(influxClientClass.instance.measurements),
                {"Project_1"},
                "Influx measurement must remain slugified during generate().",
            )



    def testExistingStateDoesNotDeleteInfluxMeasurement(self):
        existingState = (
            "[ctf]\nlastId = 8\n"
            "[gain]\nlastId = 4\n"
            "[system]\nlastId = 12\n"
        )

        with tempfile.TemporaryDirectory() as tmpDir:
            _, influxClientClass = self._createReport(
                tmpDir,
                projectName="resume-project",
                existingState=existingState,
            )

            client = influxClientClass.instance
            self.assertIsNotNone(client)
            self.assertEqual(client.switchedDatabase, "scipion")
            self.assertEqual(
                client.deletedMeasurements,
                [],
                "Continue must preserve the existing Influx measurement.",
            )

class TestInfluxTransport(unittest.TestCase):

    def testConnectPutPropagatesTransferFailure(self):
        class FailingSftp:
            def put(self, local, remote, confirm=True):
                raise OSError("simulated transfer failure")

        connect = object.__new__(Connect)
        connect.sftp = FailingSftp()
        connect.remote_path = "/remote"

        with self.assertRaisesRegex(
            OSError,
            "simulated transfer failure",
        ):
            connect.put(
                ["/local/image.jpg"],
                ["project/image.jpg"],
            )

