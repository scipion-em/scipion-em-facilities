import os
import unittest
from configparser import ConfigParser, NoSectionError
from unittest.mock import patch, MagicMock

from emfacilities import EMFACILITIES_HOME_VARNAME, Plugin
from emfacilities.protocols.report_influx import ReportInflux

class TestsInfluxReport(unittest.TestCase):
    TMP_CONF_FILE = "/tmp/kk.cfg"
    PROJ_NAME = "project"

    @patch("emfacilities.protocols.report_influx.SummaryProvider")
    def test_instantiation(self, summaryMock):

        protocol = self.getProtocolMock()

        with patch.object(ReportInflux, 'parseSecrets') as parseSecrets:
            parseSecrets.return_value = (1,2,3,4,5,6)
            with patch.object(ReportInflux, 'connectToInflux') as connectToInflux:

                # Test status config file is being generated
                report = ReportInflux(protocol, None,
                             None, None)

                parseSecrets.assert_called_once()
                connectToInflux.assert_called_once()
                self.assertTrue(os.path.exists(self.TMP_CONF_FILE), "%s tmp config file not created." % self.TMP_CONF_FILE)

                conf = ConfigParser()
                conf.read(self.TMP_CONF_FILE)
                
                self.assertEqual(conf.get("project", "projectName"),   self.PROJ_NAME, "Project name wrongly set")
                self.assertEqual(conf.get("project", "properties"), "1", "Project.properties wrongly initialized")
                self.assertEqual(conf.get("project", "acquisition"), "1", "Project.acquisition wrongly initialized")
                self.assertEqual(conf.get("project", "summary"), "1", "Project.summary wrongly initialized")
                self.assertEqual(conf.get("ctf", "lastId"), "0", "ctf.lastId wrongly initialized")
                self.assertEqual(conf.get("gain", "lastId"), "0", "gain.lastId wrongly initialized")
                self.assertEqual(conf.get("system", "lastId"), "0", "system.lastId wrongly initialized")

                # Test status config file is read
                NEW_PROJ_NAME = "newName"
                conf.set("project", "projectName", NEW_PROJ_NAME)

                with open(self.TMP_CONF_FILE, 'w') as confFile:
                    conf.write(confFile)

                # Test status config file is being read
                ReportInflux(protocol, None, None, None)

                conf = ConfigParser()
                conf.read(self.TMP_CONF_FILE)

                self.assertEqual(conf.get("project", "projectName"), NEW_PROJ_NAME, "Not reading existing status file properly")

    @classmethod
    def getProtocolMock(cls):

        protocol = MagicMock()

        if os.path.exists(cls.TMP_CONF_FILE):
            os.remove(cls.TMP_CONF_FILE)
        protocol._getTmpPath.return_value = cls.TMP_CONF_FILE
        projectMock = MagicMock()
        projectMock.getShortName.return_value = cls.PROJ_NAME
        protocol.getProject.return_value = projectMock
        return protocol

    def test_secretsParser(self):

        Plugin._defineVar(EMFACILITIES_HOME_VARNAME, "")
        protMock = self.getProtocolMock()

        with patch.object(ReportInflux, 'connectToInflux') as connectToInflux:
            # Test status config file is being read

            self.assertRaises(NoSectionError, ReportInflux, protMock, None, None, None)



if __name__ == '__main__':
    unittest.main()
