from unittest.mock import patch

from emfacilities.protocols.protocol_streaming_base import (
    ProtFacilitiesStreamingBase,
)


class _Value:
    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value


class _ClosedStreamInputSet:
    def __init__(self):
        self.uniqueCalls = []
        self.closedCalls = 0

    def getUniqueValues(self, attributes, where=None):
        self.uniqueCalls.append((attributes, where))

        if where is None:
            return list(range(1, 11))
        if where == "id > 0":
            return [9, 10]
        if where == "id > 10":
            return []

        raise AssertionError(
            "Unexpected discovery query: %r" % (where,)
        )

    def getSize(self):
        return 10

    def isStreamClosed(self):
        return True

    def close(self):
        self.closedCalls += 1


class _Image:
    def __init__(self, objId):
        self.objId = objId

    def clone(self):
        return _Image(self.objId)


class _VisibilityInputSet:
    def __init__(self, visibleIds, size):
        self._visibleIds = set(visibleIds)
        self._size = size
        self.closed = False

    def getSize(self):
        return self._size

    def __contains__(self, objId):
        return objId in self._visibleIds

    def getItem(self, field, value):
        assert field == "id"
        if value not in self._visibleIds:
            raise UnboundLocalError(
                "row not found for id %r" % value
            )
        return _Image(value)

    def close(self):
        self.closed = True


class _VisibilityOutputSet:
    STREAM_OPEN = 1
    STREAM_CLOSED = 2

    def __init__(self):
        self.ids = []

    def getSize(self):
        return len(self.ids)

    def append(self, image):
        self.ids.append(image.objId)


class _Pointer:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _ExistingOutput:
    def __init__(self):
        self.loaded = False
        self.appendEnabled = False
        self.copiedFrom = None

    def loadAllProperties(self):
        self.loaded = True

    def enableAppend(self):
        if not self.loaded:
            raise AssertionError(
                "Logical output must be refreshed before "
                "enableAppend()."
            )
        self.appendEnabled = True

    def copyInfo(self, inputs):
        self.copiedFrom = inputs

    def getSize(self):
        if not self.loaded:
            raise AssertionError(
                "Logical output must be refreshed before "
                "reading its size."
            )
        return 2

    def getIdSet(self):
        if not self.loaded:
            raise AssertionError(
                "Logical output must be refreshed before "
                "reading its ids."
            )
        return {1, 2}


class _FreshOutput:
    STREAM_OPEN = 1

    def __init__(self, filename=None):
        self.filename = filename
        self.loaded = False
        self.streamState = None
        self.copiedFrom = None

    def loadAllProperties(self):
        self.loaded = True
        raise AssertionError(
            "A backing file must not restore an output absent "
            "from protocol outputs."
        )

    def setStreamState(self, state):
        self.streamState = state

    def copyInfo(self, inputs):
        self.copiedFrom = inputs


def assert_closed_stream_reconciliation(
        testCase,
        protocolClass,
        protocolModule,
):
    class _Harness:
        def __init__(self):
            self.finished = False
            self.insertedIds = set()
            self.processedIds = set()
            self.sampleIds = set()
            self._lastInputId = 0
            self._pendingInputIds = []
            self.isStreamClosed = False
            self.lastRound = False
            self.limitReach = False
            self.timerOut = False
            self.outputSize = _Value(10)
            self.batchSize = _Value(10)
            self.inputSet = _ClosedStreamInputSet()
            self.insertedBatches = []

        def _loadInputSet(self, _):
            return self.inputSet

        # The protocol now goes through the shared discovery wrapper, so
        # the harness has to borrow it too - otherwise this test would only
        # exercise helpers nothing calls any more.
        # The stall detector lives in the shared base and runs inside
        # the reconciliation, so the harness borrows it too.
        _hasActiveStreamingWork = (
            ProtFacilitiesStreamingBase._hasActiveStreamingWork)
        _recordTerminalProgress = (
            ProtFacilitiesStreamingBase._recordTerminalProgress)

        def _discoverNewInputIds(self, knownIds):
            return ProtFacilitiesStreamingBase._discoverNewInputIds(
                self, knownIds
            )

        def _discoverIdsAfter(self, inputSet, lastId):
            return protocolClass._discoverIdsAfter(
                self,
                inputSet,
                lastId,
            )

        def _reconcileClosedStreamIds(
                self,
                inputSet,
                discoveredIds,
                knownIds,
                producerClosed,
        ):
            return protocolClass._reconcileClosedStreamIds(
                self,
                inputSet,
                discoveredIds,
                knownIds,
                producerClosed,
            )

        def isContinued(self):
            return False

        def _insertNewImageSteps(self, newIds, *args):
            ids = list(newIds)
            self.insertedBatches.append(ids)
            self.insertedIds.update(ids)
            return []

        def info(self, message):
            pass

    protocol = _Harness()

    with patch(
            protocolModule + ".time.sleep",
            return_value=None,
    ):
        protocolClass._checkNewInput(protocol)
        protocolClass._checkNewInput(protocol)

    testCase.assertEqual(
        protocol.insertedIds,
        set(range(1, 11)),
    )
    testCase.assertIn(
        ("id", None),
        protocol.inputSet.uniqueCalls,
    )


def assert_late_visibility_retry(
        testCase,
        protocolClass,
        pendingAttribute,
        processedIds,
        inputSize,
):
    class _Harness:
        def __init__(self):
            self.finished = False
            self.isStreamClosed = True
            self.processedIds = set(processedIds)
            self.sampleIds = {4, 5}
            self.limitReach = False
            self.timerOut = False
            self.outputSize = _Value(100)
            self._inputClass = object
            self._baseName = "images.sqlite"
            self.inputSet = _VisibilityInputSet(
                visibleIds={4},
                size=inputSize,
            )
            self.outputSet = _VisibilityOutputSet()
            self.streamMode = None
            self.errors = []

        def _loadInputSet(self, _):
            return self.inputSet

        def _loadOutputSet(
                self,
                SetClass,
                baseName,
                outputName=None,
        ):
            return self.outputSet

        def _updateOutputSet(
                self,
                outputName,
                outputSet,
                streamMode,
        ):
            self.streamMode = streamMode

        def _store(self):
            pass

        def error(self, msg):
            self.errors.append(msg)

    protocol = _Harness()

    protocolClass._checkNewOutput(protocol)

    testCase.assertEqual([4], protocol.outputSet.ids)
    testCase.assertEqual(
        {5},
        getattr(protocol, pendingAttribute),
    )
    testCase.assertEqual(1, len(protocol.errors))
    testCase.assertFalse(protocol.finished)
    testCase.assertEqual(
        protocol.outputSet.STREAM_OPEN,
        protocol.streamMode,
    )


def assert_persisted_output_identity(
        testCase,
        protocolClass,
        cleanPathTarget,
        outputName,
):
    inputs = object()
    existing = _ExistingOutput()

    class _Harness:
        def __init__(self):
            self.outputSet = existing
            self.inputImages = _Pointer(inputs)

        def _getPath(self, name):
            return "/tmp/" + name

    protocol = _Harness()

    loaded = protocolClass._loadOutputSet(
        protocol,
        object,
        "images.sqlite",
        outputName=outputName,
    )
    doneIds, size = protocolClass._getAllDoneIds(protocol)

    testCase.assertIs(existing, loaded)
    testCase.assertTrue(existing.loaded)
    testCase.assertTrue(existing.appendEnabled)
    testCase.assertEqual({1, 2}, set(doneIds))
    testCase.assertEqual(2, size)

    del protocol.outputSet

    with patch(cleanPathTarget) as cleanPathMock:
        fresh = protocolClass._loadOutputSet(
            protocol,
            _FreshOutput,
            "images.sqlite",
            outputName=outputName,
        )

    cleanPathMock.assert_called_once_with(
        "/tmp/images.sqlite"
    )
    testCase.assertFalse(fresh.loaded)
    testCase.assertEqual(
        _FreshOutput.STREAM_OPEN,
        fresh.streamState,
    )
    testCase.assertIs(inputs, fresh.copiedFrom)
