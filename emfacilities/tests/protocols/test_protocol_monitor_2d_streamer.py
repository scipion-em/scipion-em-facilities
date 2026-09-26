from unittest.mock import patch

from pyworkflow.tests import BaseTest, setupTestProject

from emfacilities.protocols.protocol_monitor_2d_streamer import (
    ProtMonitor2dStreamer,
)



class _Pointer:
    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value


class _InactiveInput2dProtocol:
    def __init__(self, objId=50):
        self.objId = objId

    def isActive(self):
        return False

    def getObjId(self):
        return self.objId


class _ProjectInfo:
    def getName(self):
        return "Monitor2dResumeTest"


class _OutputSet:
    def __init__(self, ids=()):
        self._ids = set(ids)

    def getIdSet(self):
        return set(self._ids)


class TestMonitor2dStreamer(BaseTest):

    def _prepareContinueMonitor(
            self,
            runIds=(),
            outputs=(),
            project=None,
    ):
        prot = self.newProtocol(
            ProtMonitor2dStreamer,
            samplingInterval=1,
        )
        prot.isContinued = lambda: True
        prot._runIds.set(list(runIds))
        prot.input2dProtocol = _Pointer(_InactiveInput2dProtocol())

        def createSubset():
            prot._counter += 1
            return object()

        prot._createSubset = createSubset
        prot.iterOutputAttributes = lambda: list(outputs)

        managerPatch = None
        if project is not None:
            prot.getProject = lambda: _ProjectInfo()
            prot.getObjId = lambda: 99

            class FakeManager:
                def loadProject(self, name):
                    return project

            managerPatch = patch(
                "emfacilities.protocols.protocol_monitor_2d_streamer.Manager",
                return_value=FakeManager(),
            )

        return prot, managerPatch


    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)


    def testContinueRestoresRunPrerequisites(self):
        prot, _ = self._prepareContinueMonitor(
            runIds=[101, 102],
        )

        restoredState = {}

        def checkNewInput():
            restoredState["prerequisites"] = list(
                prot._runPrerequisites
            )
            prot._streamClosed = True

        prot._checkNewInput = checkNewInput

        with patch(
                "emfacilities.protocols.protocol_monitor_2d_streamer.time.sleep",
                return_value=None,
        ):
            prot.monitorStep()

        self.assertEqual(
            restoredState["prerequisites"],
            [101, 102],
        )


    def testContinueRestoresProgressFromExistingSubsets(self):
        class OutputSet:
            def __init__(self, ids):
                self._ids = set(ids)

            def getIdSet(self):
                return set(self._ids)

        class Input2dProtocol:
            def isActive(self):
                return False

        class Pointer:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        prot = self.newProtocol(
            ProtMonitor2dStreamer,
            samplingInterval=1,
        )
        prot.isContinued = lambda: True
        prot.input2dProtocol = Pointer(Input2dProtocol())
        def createSubset():
            prot._counter += 1
            return object()

        prot._createSubset = createSubset
        prot.iterOutputAttributes = lambda: [
            ("outputParticles_001", OutputSet({1, 2, 3})),
            ("outputParticles_002", OutputSet({4, 5, 6})),
        ]
        prot._restoreScheduledRuns = lambda outputNames: None

        restoredState = {}

        def checkNewInput():
            restoredState["lastPartId"] = prot._lastPartId
            restoredState["counter"] = prot._counter
            restoredState["processed"] = prot._counterParticlesProcessed
            prot._streamClosed = True

        prot._checkNewInput = checkNewInput

        with patch(
                "emfacilities.protocols.protocol_monitor_2d_streamer.time.sleep",
                return_value=None,
        ):
            prot.monitorStep()

        self.assertEqual(restoredState["lastPartId"], 6)
        self.assertEqual(restoredState["counter"], 3)
        self.assertEqual(restoredState["processed"], 6)

    def testWriteSubsetPersistsScheduledRunId(self):
        class Subset:
            def write(self):
                pass

            def close(self):
                pass

        class ParticlePointer:
            def set(self, value):
                self.value = value

            def setExtended(self, value):
                self.extended = value

        class CopyProtocol:
            def __init__(self, objId):
                self.objId = objId
                self.inputParticles = ParticlePointer()

            def getObjId(self):
                return self.objId

        class Input2dProtocol:
            def getObjId(self):
                return 50

        class Pointer:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        class ProjectInfo:
            def getName(self):
                return "Monitor2dResumeTest"

        copyProt = CopyProtocol(321)

        class Project:
            def getProtocol(self, objId):
                return object()

            def copyProtocol(self, protocol):
                return copyProt

            def scheduleProtocol(self, protocol, prerequisites):
                self.prerequisites = list(prerequisites)

        project = Project()

        class FakeManager:
            def loadProject(self, name):
                return project

        prot = self.newProtocol(ProtMonitor2dStreamer)
        prot._counter = 3
        prot._runPrerequisites = [101, 102]
        prot.input2dProtocol = Pointer(Input2dProtocol())
        prot.getProject = lambda: ProjectInfo()
        prot.getObjId = lambda: 99
        prot._defineOutputs = lambda **kwargs: None
        prot._defineTransformRelation = lambda *args: None

        stored = []
        prot._store = lambda obj=None: stored.append(obj)

        with patch(
                "emfacilities.protocols.protocol_monitor_2d_streamer.Manager",
                return_value=FakeManager(),
        ):
            prot._writeSubset(Subset())

        self.assertEqual(list(prot._runIds), [321])
        self.assertTrue(any(obj is prot._runIds for obj in stored))
        self.assertEqual(prot._runPrerequisites, [101, 102, 321])


    def testContinueRecoversScheduledRunMissingFromRunIds(self):
        class ParentProtocol:
            def getObjId(self):
                return 99

        class ParticlePointer:
            def getObjValue(self):
                return ParentProtocol()

            def getExtended(self):
                return "outputParticles_001"

        class ScheduledRun:
            inputParticles = ParticlePointer()

            def getObjId(self):
                return 321

        class Project:
            def getRuns(self):
                return [ScheduledRun()]

        prot, managerPatch = self._prepareContinueMonitor(
            outputs=[
                ("outputParticles_001", _OutputSet({1, 2})),
            ],
            project=Project(),
        )

        restoredState = {}

        def checkNewInput():
            restoredState["runIds"] = list(prot._runIds)
            restoredState["prerequisites"] = list(
                prot._runPrerequisites
            )
            prot._streamClosed = True

        prot._checkNewInput = checkNewInput

        with managerPatch:
            prot.monitorStep()

        self.assertEqual(restoredState["runIds"], [321])
        self.assertEqual(
            restoredState["prerequisites"],
            [321],
        )



    def testContinueSchedulesRunForPersistedSubsetMissingClassification(self):
        class ParentProtocol:
            def getObjId(self):
                return 99

        class ParticlePointer:
            def __init__(self):
                self.parent = None
                self.extended = None

            def set(self, value):
                self.parent = value

            def setExtended(self, value):
                self.extended = value

        class CopyProtocol:
            def __init__(self):
                self.inputParticles = ParticlePointer()

            def getObjId(self):
                return 321

        class Project:
            def __init__(self):
                self.scheduled = []

            def getRuns(self):
                return []

            def getProtocol(self, objId):
                if objId == 99:
                    return ParentProtocol()
                return object()

            def copyProtocol(self, protocol):
                return CopyProtocol()

            def scheduleProtocol(self, protocol, prerequisites):
                self.scheduled.append(
                    (
                        protocol.inputParticles.extended,
                        list(prerequisites),
                    )
                )

        project = Project()
        prot, managerPatch = self._prepareContinueMonitor(
            outputs=[
                ("outputParticles_001", _OutputSet({1, 2})),
            ],
            project=project,
        )

        restoredState = {}

        def checkNewInput():
            restoredState["runIds"] = list(prot._runIds)
            restoredState["prerequisites"] = list(
                prot._runPrerequisites
            )
            prot._streamClosed = True

        prot._checkNewInput = checkNewInput

        with managerPatch:
            prot.monitorStep()

        self.assertEqual(
            project.scheduled,
            [("outputParticles_001", [])],
        )
        self.assertEqual(restoredState["runIds"], [321])
        self.assertEqual(
            restoredState["prerequisites"],
            [321],
        )


    def testStartingNumberSkipsParticleCountNotIds(self):
        class Particle:
            def __init__(self, objId):
                self._objId = objId

            def getObjId(self):
                return self._objId

        class InputSet:
            def load(self):
                pass

            def loadAllProperties(self):
                pass

            def isStreamClosed(self):
                return False

            def iterItems(self, **kwargs):
                return iter([
                    Particle(10),
                    Particle(20),
                    Particle(30),
                    Particle(40),
                ])

            def close(self):
                pass

        prot = self.newProtocol(
            ProtMonitor2dStreamer,
            startingNumber=2,
        )
        prot.inputParticles = _Pointer(InputSet())
        prot._lastPartId = 0
        prot._streamClosed = False

        particleIds = [
            particle.getObjId()
            for particle in prot._iterParticles()
        ]

        self.assertEqual(particleIds, [30, 40])

    def testMonitorDoesNotSleepAfterStreamCloses(self):
        prot = self.newProtocol(
            ProtMonitor2dStreamer,
            samplingInterval=10,
        )
        prot.isContinued = lambda: False
        prot.input2dProtocol = _Pointer(
            _InactiveInput2dProtocol()
        )
        prot._createSubset = lambda: object()

        def checkNewInput():
            prot._streamClosed = True

        prot._checkNewInput = checkNewInput

        with patch(
                "emfacilities.protocols.protocol_monitor_2d_streamer.time.sleep",
        ) as sleepMock:
            prot.monitorStep()

        sleepMock.assert_not_called()

    def testClosedStreamWithNoNewParticlesDoesNotWriteEmptySubset(self):
        class EmptySubset:
            def getSize(self):
                return 0

        prot = self.newProtocol(ProtMonitor2dStreamer)
        prot._subset = EmptySubset()
        prot._lastMicId = None
        prot._lastPartId = 6
        prot._counterNewParticles = 0
        prot._counterParticlesProcessed = 6
        prot._streamClosed = False

        def iterParticles():
            prot._streamClosed = True
            return iter(())

        prot._iterParticles = iterParticles

        written = []
        prot._writeSubset = lambda subset: written.append(subset)

        prot._checkNewInput()

        self.assertEqual(written, [])

    def testParticleLimitStopsAtExactMaximum(self):
        prot = self.newProtocol(
            ProtMonitor2dStreamer,
            maximumOption=ProtMonitor2dStreamer.NUMBER_PARTICLES,
            numberParticles=100,
        )
        prot._counterParticlesProcessed = 100

        self.assertTrue(prot.classificationStop())

    def testParticleLimitWritesPendingSubsetWithoutExtraParticle(self):
        class Particle:
            def __init__(self, objId, micId):
                self._objId = objId
                self._micId = micId

            def getObjId(self):
                return self._objId

            def getMicId(self):
                return self._micId

        class Subset:
            def __init__(self):
                self.ids = []

            def append(self, particle):
                self.ids.append(particle.getObjId())

            def getSize(self):
                return len(self.ids)

        prot = self.newProtocol(
            ProtMonitor2dStreamer,
            maximumOption=ProtMonitor2dStreamer.NUMBER_PARTICLES,
            numberParticles=2,
            batchSize=100,
        )
        prot._subset = Subset()
        prot._lastMicId = None
        prot._lastPartId = 0
        prot._counterNewParticles = 0
        prot._counterParticlesProcessed = 0
        prot._streamClosed = False
        prot._iterParticles = lambda: iter([
            Particle(1, 1),
            Particle(2, 1),
            Particle(3, 2),
        ])

        written = []
        prot._writeSubset = lambda subset: written.append(list(subset.ids))

        prot._checkNewInput()

        self.assertTrue(prot._streamClosed)
        self.assertEqual(written, [[1, 2]])
        self.assertEqual(prot._lastPartId, 2)
        self.assertEqual(prot._counterParticlesProcessed, 2)

    def testBatchBoundaryDoesNotSplitNextMicrograph(self):
        class Particle:
            def __init__(self, objId, micId):
                self._objId = objId
                self._micId = micId

            def getObjId(self):
                return self._objId

            def getMicId(self):
                return self._micId

        class Subset:
            def __init__(self):
                self.ids = []

            def append(self, particle):
                self.ids.append(particle.getObjId())

            def getSize(self):
                return len(self.ids)

        prot = self.newProtocol(
            ProtMonitor2dStreamer,
            batchSize=2,
        )
        prot._subset = Subset()
        prot._lastMicId = None
        prot._lastPartId = 0
        prot._counterNewParticles = 0
        prot._counterParticlesProcessed = 0
        prot._streamClosed = False
        prot.maximumOption.set(ProtMonitor2dStreamer.NONE_OPTION)

        particles = [
            Particle(1, 1),
            Particle(2, 1),
            Particle(3, 1),
            Particle(4, 2),
            Particle(5, 2),
        ]

        prot._iterParticles = lambda: iter(particles)

        written = []

        def writeSubset(subset):
            written.append(list(subset.ids))

        def createSubset():
            return Subset()

        prot._writeSubset = writeSubset
        prot._createSubset = createSubset

        prot._checkNewInput()

        self.assertEqual(
            written,
            [[1, 2, 3]],
            "The first particle from the next micrograph must not be "
            "written into the previous batch.",
        )
        self.assertEqual(
            prot._subset.ids,
            [4, 5],
            "The next micrograph must remain intact in the new subset.",
        )
