import unittest
from unittest.mock import patch

from emfacilities.protocols.report_html import (
    MIC_ID,
    MIC_PATH,
    MIC_THUMBS,
    PSD_PATH,
    PSD_THUMBS,
    SHIFT_PATH,
    SHIFT_THUMBS,
    ReportHtml,
)


class TestReportHtml(unittest.TestCase):

    def testCheckNewThumbsReadyStopsAtFirstMissingThumbnail(self):
        import os
        import tempfile

        with tempfile.TemporaryDirectory() as tmpDir:
            report = object.__new__(ReportHtml)
            report.reportDir = tmpDir
            report.thumbsReady = 0
            report.thumbPaths = {
                MIC_THUMBS: [
                    "imgMicThumbs/mic001.jpg",
                    "imgMicThumbs/mic002.jpg",
                    "imgMicThumbs/mic003.jpg",
                ],
                PSD_THUMBS: [
                    "imgPsdThumbs/psd001.jpg",
                    "imgPsdThumbs/psd002.jpg",
                    "imgPsdThumbs/psd003.jpg",
                ],
            }

            for relPath in (
                "imgMicThumbs/mic002.jpg",
                "imgMicThumbs/mic003.jpg",
                "imgPsdThumbs/psd002.jpg",
                "imgPsdThumbs/psd003.jpg",
            ):
                absPath = os.path.join(tmpDir, relPath)
                os.makedirs(os.path.dirname(absPath), exist_ok=True)
                open(absPath, "w").close()

            self.assertEqual(
                report.checkNewThumbsReady(),
                0,
                "Ready thumbnails after a missing earlier thumbnail must "
                "not advance the contiguous ready prefix.",
            )

            for relPath in (
                "imgMicThumbs/mic001.jpg",
                "imgPsdThumbs/psd001.jpg",
            ):
                absPath = os.path.join(tmpDir, relPath)
                open(absPath, "w").close()

            self.assertEqual(
                report.checkNewThumbsReady(),
                3,
            )




    def testGetThumbPathsDropsPsdKeysWhenMicrographHasNoPsd(self):
        class Micrograph:
            def getFileName(self):
                return "/tmp/mic002.mrc"

        class OutputSet:
            def getIdSet(self):
                return {2}

            def __getitem__(self, itemId):
                return Micrograph()

        class AlignProtocol:
            outputMicrographs = OutputSet()

        report = object.__new__(ReportHtml)
        report.alignProtocol = AlignProtocol()
        report.ctfProtocol = None
        report.thumbPaths = {
            MIC_THUMBS: [],
            PSD_THUMBS: [],
            SHIFT_THUMBS: [],
            MIC_PATH: [],
            SHIFT_PATH: [],
            PSD_PATH: [],
            MIC_ID: [],
        }

        with patch(
                "emfacilities.protocols.report_html.getUpdatedProtocol",
                side_effect=lambda protocol: protocol,
        ):
            report.getThumbPaths()

        self.assertNotIn(
            PSD_THUMBS,
            report.thumbPaths,
            "Missing PSD data must not create a fake PSD thumbnail.",
        )
        self.assertNotIn(PSD_PATH, report.thumbPaths)
        self.assertNotIn(
            "None.jpg",
            str(report.thumbPaths),
        )


    def testGetThumbPathsKeepsExistingPsdWithoutCtfProtocol(self):
        class FileRef:
            def __init__(self, fileName):
                self.fileName = fileName

            def getFileName(self):
                return self.fileName

        class Micrograph:
            psdJpeg = FileRef("/tmp/mic001_psd.jpg")

            def getFileName(self):
                return "/tmp/mic001.mrc"

        class OutputSet:
            def getIdSet(self):
                return {1}

            def __getitem__(self, itemId):
                self.lastItemId = itemId
                return Micrograph()

        class AlignProtocol:
            outputMicrographs = OutputSet()

        report = object.__new__(ReportHtml)
        report.alignProtocol = AlignProtocol()
        report.ctfProtocol = None
        report.thumbPaths = {
            MIC_THUMBS: [],
            PSD_THUMBS: [],
            SHIFT_THUMBS: [],
            MIC_PATH: [],
            SHIFT_PATH: [],
            PSD_PATH: [],
            MIC_ID: [],
        }

        with patch(
                "emfacilities.protocols.report_html.getUpdatedProtocol",
                side_effect=lambda protocol: protocol,
        ):
            report.getThumbPaths()

        self.assertIn(
            PSD_THUMBS,
            report.thumbPaths,
            "An existing micrograph PSD must remain available to the report.",
        )
        self.assertIn(PSD_PATH, report.thumbPaths)
        self.assertEqual(
            report.thumbPaths[PSD_PATH],
            ["/tmp/mic001_psd.jpg"],
        )
        self.assertEqual(
            report.thumbPaths[PSD_THUMBS],
            ["imgPsdThumbs/mic001_psd.jpg"],
        )
