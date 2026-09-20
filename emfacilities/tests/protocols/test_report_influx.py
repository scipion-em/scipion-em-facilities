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


class TestReportInfluxResume(unittest.TestCase):

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
        class Project:
            def getShortName(self):
                return "transfer-project"

            def getCreationTime(self):
                return datetime(2026, 1, 1, 12, 0, 0)

        class Protocol:
            def __init__(self, workingDir):
                self.workingDir = workingDir
                self.reportPath = os.path.join(workingDir, "index.html")
                self.reportDir = workingDir
                self.project = Project()

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

        class Provider:
            acquisition = []

            def refreshObjects(self):
                pass

            def getObjects(self):
                return []

        class SystemMonitor:
            def getData(self, lastId=-1):
                return [
                    {
                        "id": 1,
                        "timestamp": datetime(2026, 1, 1, 12, 0, 1),
                        "cpu": 10.0,
                    }
                ]

        class FakeInfluxDBClient:
            def __init__(self, **kwargs):
                pass

            def switch_database(self, name):
                pass

            def delete_series(self, measurement):
                pass

            def create_retention_policy(self, *args, **kwargs):
                pass

            def write_points(self, points):
                pass

        with tempfile.TemporaryDirectory() as tmpDir:
            with open(os.path.join(tmpDir, "secrets.cfg"), "w") as handle:
                handle.write(
                    "[influx]\n"
                    "dataBase = scipion\n"
                    "passwordInflux = cGFzcw==\n"
                    "usernameInflux = dXNlcg==\n"
                    "hostinflux = localhost\n"
                    "port = 8086\n"
                    "ssl = false\n"
                    "verify_ssl = false\n"
                    "TimeDelta = 0\n"
                    "apacheImgDir = /tmp\n"
                    "[paramiko]\n"
                    "usernameParamiko = user\n"
                    "passwordParamiko = pass\n"
                    "keyfilepath = key\n"
                    "keyfiletype = rsa\n"
                    "remote_path = /tmp\n"
                    "hostparamiko = localhost\n"
                )

            fakeInfluxModule = types.ModuleType("influxdb")
            fakeInfluxModule.InfluxDBClient = FakeInfluxDBClient

            with patch.dict(
                os.environ,
                {
                    EMFACILITIES_HOME_VARNAME: tmpDir,
                    "SCIPION_VERSION": "test",
                },
            ), patch.dict(
                sys.modules,
                {"influxdb": fakeInfluxModule},
            ):
                report = ReportInflux(
                    Protocol(tmpDir),
                    ctfMonitor=None,
                    sysMonitor=SystemMonitor(),
                    movieGainMonitor=None,
                )
                report.provider = Provider()
                report.transferFiles = lambda: False

                self.assertFalse(
                    report.generate(finished=False),
                    "A failed image transfer must keep the summary monitor running.",
                )



    def testGainProgressIsPersistedAfterEachSuccessfulWrite(self):
        class Project:
            def getShortName(self):
                return "gain-resume-project"

            def getCreationTime(self):
                return datetime(2026, 1, 1, 12, 0, 0)

        class Protocol:
            def __init__(self, workingDir):
                self.workingDir = workingDir
                self.reportPath = os.path.join(workingDir, "index.html")
                self.reportDir = workingDir
                self.project = Project()

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

        class Provider:
            acquisition = []

            def refreshObjects(self):
                pass

            def getObjects(self):
                return []

        class MovieGainMonitor:
            def getData(self, lastId=-1):
                return [
                    {"idx": 1, "stddev": 0.1, "ratio1": 1.1, "ratio2": 1.2},
                    {"idx": 2, "stddev": 0.2, "ratio1": 1.3, "ratio2": 1.4},
                ]

        class FakeInfluxDBClient:
            instance = None

            def __init__(self, **kwargs):
                type(self).instance = self
                self.gainWrites = 0

            def switch_database(self, name):
                pass

            def delete_series(self, measurement):
                pass

            def create_retention_policy(self, *args, **kwargs):
                pass

            def write_points(self, points):
                point = points[0]
                if point.get("tags", {}).get("section") == "gain":
                    self.gainWrites += 1
                    if self.gainWrites == 2:
                        raise RuntimeError("simulated influx failure")

        with tempfile.TemporaryDirectory() as tmpDir:
            with open(os.path.join(tmpDir, "secrets.cfg"), "w") as handle:
                handle.write(
                    "[influx]\n"
                    "dataBase = scipion\n"
                    "passwordInflux = cGFzcw==\n"
                    "usernameInflux = dXNlcg==\n"
                    "hostinflux = localhost\n"
                    "port = 8086\n"
                    "ssl = false\n"
                    "verify_ssl = false\n"
                    "TimeDelta = 0\n"
                    "apacheImgDir = /tmp\n"
                    "[paramiko]\n"
                    "usernameParamiko = user\n"
                    "passwordParamiko = pass\n"
                    "keyfilepath = key\n"
                    "keyfiletype = rsa\n"
                    "remote_path = /tmp\n"
                    "hostparamiko = localhost\n"
                )

            fakeInfluxModule = types.ModuleType("influxdb")
            fakeInfluxModule.InfluxDBClient = FakeInfluxDBClient

            with patch.dict(
                os.environ,
                {
                    EMFACILITIES_HOME_VARNAME: tmpDir,
                    "SCIPION_VERSION": "test",
                },
            ), patch.dict(
                sys.modules,
                {"influxdb": fakeInfluxModule},
            ):
                report = ReportInflux(
                    Protocol(tmpDir),
                    ctfMonitor=None,
                    sysMonitor=None,
                    movieGainMonitor=MovieGainMonitor(),
                )
                report.provider = Provider()
                report.transferFiles = lambda: True

                with self.assertRaisesRegex(
                    RuntimeError,
                    "simulated influx failure",
                ):
                    report.generate(finished=False)

            report.confParser.read(report.confFileName)
            self.assertEqual(
                report.confParser.getint("gain", "lastId"),
                1,
                "A successful gain write must be checkpointed before the next item.",
            )



    def testGainPointsUseDistinctTimestampsWithinSameBatch(self):
        class Project:
            def getShortName(self):
                return "gain-project"

            def getCreationTime(self):
                return datetime(2026, 1, 1, 12, 0, 0)

        class Protocol:
            def __init__(self, workingDir):
                self.workingDir = workingDir
                self.reportPath = os.path.join(workingDir, "index.html")
                self.reportDir = workingDir
                self.project = Project()

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

        class Provider:
            acquisition = []

            def refreshObjects(self):
                pass

            def getObjects(self):
                return []

        class MovieGainMonitor:
            def getData(self, lastId=-1):
                return [
                    {"idx": 1, "stddev": 0.1, "ratio1": 1.1, "ratio2": 1.2},
                    {"idx": 2, "stddev": 0.2, "ratio1": 1.3, "ratio2": 1.4},
                ]

        class FakeInfluxDBClient:
            instance = None

            def __init__(self, **kwargs):
                type(self).instance = self
                self.gainTimes = []

            def switch_database(self, name):
                pass

            def delete_series(self, measurement):
                pass

            def create_retention_policy(self, *args, **kwargs):
                pass

            def write_points(self, points):
                for point in points:
                    if point.get("tags", {}).get("section") == "gain":
                        self.gainTimes.append(point["time"])

        with tempfile.TemporaryDirectory() as tmpDir:
            with open(os.path.join(tmpDir, "secrets.cfg"), "w") as handle:
                handle.write(
                    "[influx]\n"
                    "dataBase = scipion\n"
                    "passwordInflux = cGFzcw==\n"
                    "usernameInflux = dXNlcg==\n"
                    "hostinflux = localhost\n"
                    "port = 8086\n"
                    "ssl = false\n"
                    "verify_ssl = false\n"
                    "TimeDelta = 0\n"
                    "apacheImgDir = /tmp\n"
                    "[paramiko]\n"
                    "usernameParamiko = user\n"
                    "passwordParamiko = pass\n"
                    "keyfilepath = key\n"
                    "keyfiletype = rsa\n"
                    "remote_path = /tmp\n"
                    "hostparamiko = localhost\n"
                )

            fakeInfluxModule = types.ModuleType("influxdb")
            fakeInfluxModule.InfluxDBClient = FakeInfluxDBClient

            with patch.dict(
                os.environ,
                {
                    EMFACILITIES_HOME_VARNAME: tmpDir,
                    "SCIPION_VERSION": "test",
                },
            ), patch.dict(
                sys.modules,
                {"influxdb": fakeInfluxModule},
            ):
                report = ReportInflux(
                    Protocol(tmpDir),
                    ctfMonitor=None,
                    sysMonitor=None,
                    movieGainMonitor=MovieGainMonitor(),
                )
                report.provider = Provider()
                report.transferFiles = lambda: True
                report.generate(finished=False)

            gainTimes = FakeInfluxDBClient.instance.gainTimes
            self.assertEqual(len(gainTimes), 2)
            self.assertNotEqual(
                gainTimes[0],
                gainTimes[1],
                "Gain points in the same batch must not share an Influx timestamp.",
            )



    def testGenerateKeepsSlugifiedMeasurementName(self):
        class Project:
            def getShortName(self):
                return "Project #1"

            def getCreationTime(self):
                return datetime(2026, 1, 1, 12, 0, 0)

        class Protocol:
            def __init__(self, workingDir):
                self.workingDir = workingDir
                self.reportPath = os.path.join(workingDir, "index.html")
                self.reportDir = workingDir
                self.project = Project()

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

        class Provider:
            acquisition = []

            def refreshObjects(self):
                pass

            def getObjects(self):
                return []

        class FakeInfluxDBClient:
            instance = None

            def __init__(self, **kwargs):
                type(self).instance = self
                self.measurements = []

            def switch_database(self, name):
                pass

            def delete_series(self, measurement):
                pass

            def create_retention_policy(self, *args, **kwargs):
                pass

            def write_points(self, points):
                self.measurements.extend(
                    point["measurement"] for point in points
                )

        with tempfile.TemporaryDirectory() as tmpDir:
            with open(os.path.join(tmpDir, "secrets.cfg"), "w") as handle:
                handle.write(
                    "[influx]\n"
                    "dataBase = scipion\n"
                    "passwordInflux = cGFzcw==\n"
                    "usernameInflux = dXNlcg==\n"
                    "hostinflux = localhost\n"
                    "port = 8086\n"
                    "ssl = false\n"
                    "verify_ssl = false\n"
                    "TimeDelta = 0\n"
                    "apacheImgDir = /tmp\n"
                    "[paramiko]\n"
                    "usernameParamiko = user\n"
                    "passwordParamiko = pass\n"
                    "keyfilepath = key\n"
                    "keyfiletype = rsa\n"
                    "remote_path = /tmp\n"
                    "hostparamiko = localhost\n"
                )

            fakeInfluxModule = types.ModuleType("influxdb")
            fakeInfluxModule.InfluxDBClient = FakeInfluxDBClient

            with patch.dict(
                os.environ,
                {
                    EMFACILITIES_HOME_VARNAME: tmpDir,
                    "SCIPION_VERSION": "test",
                },
            ), patch.dict(
                sys.modules,
                {"influxdb": fakeInfluxModule},
            ):
                report = ReportInflux(
                    Protocol(tmpDir),
                    ctfMonitor=None,
                    sysMonitor=None,
                    movieGainMonitor=None,
                )
                report.provider = Provider()
                report.transferFiles = lambda: True
                report.generate(finished=False)

            self.assertTrue(FakeInfluxDBClient.instance.measurements)
            self.assertEqual(
                set(FakeInfluxDBClient.instance.measurements),
                {"Project_1"},
                "Influx measurement must remain slugified during generate().",
            )



    def testExistingStateDoesNotDeleteInfluxMeasurement(self):
        class Project:
            def getShortName(self):
                return "resume-project"

        class Protocol:
            def __init__(self, workingDir):
                self.workingDir = workingDir
                self.reportPath = os.path.join(workingDir, "index.html")
                self.reportDir = workingDir

            def _getCtfProtocol(self):
                return None

            def _getAlignProtocol(self):
                return None

            def _getTmpPath(self, name):
                return os.path.join(self.workingDir, name)

            def getProject(self):
                return Project()

            def getInputProtocols(self):
                return []

        class FakeInfluxDBClient:
            instance = None

            def __init__(self, **kwargs):
                type(self).instance = self
                self.switchedDatabase = None
                self.deletedMeasurements = []

            def switch_database(self, name):
                self.switchedDatabase = name

            def delete_series(self, measurement):
                self.deletedMeasurements.append(measurement)

            def create_retention_policy(self, *args, **kwargs):
                pass

        with tempfile.TemporaryDirectory() as tmpDir:
            with open(os.path.join(tmpDir, CONFILE), "w") as handle:
                handle.write(
                    "[ctf]\nlastId = 8\n"
                    "[gain]\nlastId = 4\n"
                    "[system]\nlastId = 12\n"
                )

            with open(os.path.join(tmpDir, "secrets.cfg"), "w") as handle:
                handle.write(
                    "[influx]\n"
                    "dataBase = scipion\n"
                    "passwordInflux = cGFzcw==\n"
                    "usernameInflux = dXNlcg==\n"
                    "hostinflux = localhost\n"
                    "port = 8086\n"
                    "ssl = false\n"
                    "verify_ssl = false\n"
                    "TimeDelta = 0\n"
                    "apacheImgDir = /tmp\n"
                    "[paramiko]\n"
                    "usernameParamiko = user\n"
                    "passwordParamiko = pass\n"
                    "keyfilepath = key\n"
                    "keyfiletype = rsa\n"
                    "remote_path = /tmp\n"
                    "hostparamiko = localhost\n"
                )

            fakeInfluxModule = types.ModuleType("influxdb")
            fakeInfluxModule.InfluxDBClient = FakeInfluxDBClient

            with patch.dict(
                os.environ,
                {EMFACILITIES_HOME_VARNAME: tmpDir},
            ), patch.dict(
                sys.modules,
                {"influxdb": fakeInfluxModule},
            ):
                ReportInflux(
                    Protocol(tmpDir),
                    ctfMonitor=None,
                    sysMonitor=None,
                    movieGainMonitor=None,
                )

            client = FakeInfluxDBClient.instance
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

