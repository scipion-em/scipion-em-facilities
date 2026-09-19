import os
import sys
import tempfile
import types
import unittest
from datetime import datetime
from unittest.mock import patch

from emfacilities.constants import EMFACILITIES_HOME_VARNAME
from emfacilities.protocols.report_influx import CONFILE, ReportInflux


class TestReportInfluxResume(unittest.TestCase):

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
