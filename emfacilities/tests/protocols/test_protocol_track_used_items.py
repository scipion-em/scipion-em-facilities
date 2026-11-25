import unittest

from emfacilities.protocols.protocol_trackUsedItems import UsedItemsTracker


class TestUsedItemsTrackerStepSelection(unittest.TestCase):

    def testSaveSetAsJPGRetryToleratesExistingOutputDirectory(self):
        import os
        import tempfile

        tracker = object.__new__(UsedItemsTracker)

        with tempfile.TemporaryDirectory() as tmpDir:
            outDir = os.path.join(tmpDir, "usedParticles")
            os.mkdir(outDir)

            tracker.saveSetAsJPG([], outDir)

        self.assertTrue(True)



    def testNoiseCoordinatesUseDimensionsFromMatchingMicrograph(self):
        import os
        import tempfile
        from unittest.mock import patch

        tracker = object.__new__(UsedItemsTracker)
        object.__setattr__(
            tracker,
            "micDic",
            {
                1: "/tmp/mic001.mrc",
                2: "/tmp/mic002.mrc",
            },
        )
        object.__setattr__(tracker, "boxSize", 64)

        class IntParam:
            def get(self):
                return 1

        object.__setattr__(tracker, "numberOfThreads", IntParam())

        dimensions = {
            "/tmp/mic001.mrc": (100, 110, 1, 1),
            "/tmp/mic002.mrc": (200, 210, 1, 1),
        }

        class FakeImage:
            def __init__(self, dims):
                self.dims = dims

            def getDimensions(self):
                return self.dims

        class FakeImageHandler:
            def read(self, fileName):
                return FakeImage(dimensions[fileName])

        captured = []

        class FakeParallel:
            def __init__(self, **kwargs):
                pass

            def __call__(self, jobs):
                captured.extend(list(jobs))

        def fakeDelayed(function):
            def buildJob(*args):
                return args
            return buildJob

        with tempfile.TemporaryDirectory() as coordsDir:
            open(os.path.join(coordsDir, "mic001.pos"), "w").close()
            open(os.path.join(coordsDir, "mic002.pos"), "w").close()

            with patch(
                "emfacilities.protocols.protocol_trackUsedItems.ih",
                return_value=FakeImageHandler(),
            ), patch(
                "emfacilities.protocols.protocol_trackUsedItems.Parallel",
                FakeParallel,
            ), patch(
                "emfacilities.protocols.protocol_trackUsedItems.delayed",
                fakeDelayed,
            ):
                tracker.getNoiseCoordinates(
                    coordsDir,
                    "/tmp/noise",
                    noiseNumber=5,
                )

        dimsByMic = {
            job[1]: job[3]
            for job in captured
        }
        self.assertEqual(
            dimsByMic["/tmp/mic001.mrc"],
            (100, 110),
        )
        self.assertEqual(
            dimsByMic["/tmp/mic002.mrc"],
            (200, 210),
        )



    def testOutputGraphIsNotSharedAcrossCalls(self):
        first = object.__new__(UsedItemsTracker)
        object.__setattr__(first, "inpDic", {1: {}})
        firstGraph = first.generateOutputsGraphRec(
            1,
            prevCode="first.output",
        )

        second = object.__new__(UsedItemsTracker)
        object.__setattr__(second, "inpDic", {2: {}})
        secondGraph = second.generateOutputsGraphRec(
            2,
            prevCode="second.output",
        )

        self.assertIn("first.output", firstGraph)
        self.assertIn("second.output", secondGraph)
        self.assertNotIn(
            "first.output",
            secondGraph,
            "Each graph generation must start from a fresh graph.",
        )



    def testClassTrackingStepsFollowTheirOwnFlags(self):
        class BoolParam:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        def scheduledSteps(track2D, track3D):
            tracker = object.__new__(UsedItemsTracker)
            object.__setattr__(tracker, "trackParticles", BoolParam(False))
            object.__setattr__(tracker, "trackMics", BoolParam(False))
            object.__setattr__(tracker, "trackClasses2D", BoolParam(track2D))
            object.__setattr__(tracker, "trackClasses3D", BoolParam(track3D))
            object.__setattr__(tracker, "saveJPG", BoolParam(False))

            calls = []

            def insertStep(name, prerequisites=None):
                calls.append(name)
                return len(calls)

            object.__setattr__(tracker, "_insertFunctionStep", insertStep)
            tracker._insertAllSteps()
            return calls

        only3D = scheduledSteps(track2D=False, track3D=True)
        self.assertIn("trackClasses3DStep", only3D)
        self.assertNotIn("trackClasses2DStep", only3D)

        only2D = scheduledSteps(track2D=True, track3D=False)
        self.assertIn("trackClasses2DStep", only2D)
        self.assertNotIn("trackClasses3DStep", only2D)
