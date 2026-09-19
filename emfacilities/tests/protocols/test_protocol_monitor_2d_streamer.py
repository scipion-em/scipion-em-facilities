from unittest.mock import patch

from pyworkflow.tests import BaseTest, setupTestProject

from emfacilities.protocols.protocol_monitor_2d_streamer import (
    ProtMonitor2dStreamer,
)


class TestMonitor2dStreamer(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)


    def testContinueRestoresRunPrerequisites(self):
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
        prot._runIds.set([101, 102])
        prot.input2dProtocol = Pointer(Input2dProtocol())
        prot._createSubset = lambda: object()
        prot.iterOutputAttributes = lambda: []

        restoredState = {}

        def checkNewInput():
            restoredState["prerequisites"] = list(prot._runPrerequisites)
            prot._streamClosed = True

        prot._checkNewInput = checkNewInput

        with patch(
                "emfacilities.protocols.protocol_monitor_2d_streamer.time.sleep",
                return_value=None,
        ):
            prot.monitorStep()

        self.assertEqual(restoredState["prerequisites"], [101, 102])


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

        class Pointer:
            def __init__(self, value):
                self.value = value

            def get(self):
                return self.value

        prot = self.newProtocol(
            ProtMonitor2dStreamer,
            startingNumber=2,
        )
        prot.inputParticles = Pointer(InputSet())
        prot._lastPartId = 0
        prot._streamClosed = False

        particleIds = [p.getObjId() for p in prot._iterParticles()]

        self.assertEqual(particleIds, [30, 40])



