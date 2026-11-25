# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import pyworkflow.protocol.constants as cons
from pyworkflow.protocol import ProtStreamingBase
import pyworkflow.utils as pwutils
from pwem.protocols import EMProtocol


class ProtFacilitiesStreamingBase(EMProtocol, ProtStreamingBase):
    _label = None

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    def _loadLogicalSet(self, pointer):
        inputSet = pointer.get()
        inputSet.loadAllProperties()
        return inputSet

    def _discoverIdsAfter(self, inputSet, lastId):
        ids = list(inputSet.getUniqueValues('id', where='id > %d' % lastId))
        if ids:
            lastId = max(ids)
        return ids, lastId

    # How many polls a closed producer may keep showing exactly the same
    # incomplete view before the protocol gives up on it.
    TERMINAL_STALL_POLLS = 10

    def _hasActiveStreamingWork(self):
        """Whether something is still in flight for this protocol.

        Work in flight is progress, however long it takes, so it must
        never be counted towards a terminal stall. A protocol that does
        not know how to answer this says so by returning True, which
        simply means the stall detector stays out of its way.
        """
        return True

    def _recordTerminalProgress(self, inputSet, knownIds, watermarkAttr,
                                terminalConsistent):
        """Refuse to poll forever for rows that are never coming.

        A producer can close declaring more items than the consumer can
        see, and usually the rest turn up a moment later. When they do
        not - the declared size, what is known, and the watermark all
        stay exactly as they were, poll after poll, with nothing in
        flight - the protocol would otherwise sit there RUNNING for the
        rest of time. Say what is missing and fail instead.
        """
        if terminalConsistent:
            self._terminalStallSignature = None
            self._terminalStallCount = 0

            return

        if self._hasActiveStreamingWork():
            self._terminalStallCount = 0

            return

        signature = (inputSet.getSize(), len(knownIds),
                     getattr(self, watermarkAttr, 0))

        if signature == getattr(self, '_terminalStallSignature', None):
            self._terminalStallCount = getattr(
                self, '_terminalStallCount', 0) + 1
        else:
            self._terminalStallSignature = signature
            self._terminalStallCount = 1

        if self._terminalStallCount >= self.TERMINAL_STALL_POLLS:
            raise RuntimeError(
                "The input stream closed declaring %d items but only %d "
                "are visible, and that has not changed in %d polls with "
                "nothing left to process. Refusing to wait for rows that "
                "are not coming."
                % (inputSet.getSize(), len(knownIds),
                   self._terminalStallCount))

    def _reconcileClosedStreamIds(
            self,
            inputSet,
            discoveredIds,
            knownIds,
            producerClosed,
    ):
        """Reconcile IDs once a producer closes but rows lag behind metadata.

        Normal streaming discovery remains monotonic (``id > watermark``).
        A Set can transiently expose a closed stream whose declared size is
        ahead of the rows visible to the consumer. In that terminal mismatch
        only, scan the logical IDs again so late-visible rows below the
        watermark are not lost forever.

        Returns ``(newIds, terminalConsistent)``. ``terminalConsistent`` is
        False only while a closed producer still declares more items than are
        currently visible to this consumer.
        """
        discoveredIds = list(discoveredIds)

        if not producerClosed:
            return discoveredIds, True

        expectedSize = inputSet.getSize()
        knownIds = set(knownIds)
        visibleKnownIds = knownIds.union(discoveredIds)

        if len(visibleKnownIds) >= expectedSize:
            self._recordTerminalProgress(inputSet, visibleKnownIds,
                                         '_lastInputId', True)

            return discoveredIds, True

        reconciledIds = list(
            inputSet.getUniqueValues('id')
        )

        if reconciledIds:
            self._lastInputId = max(
                self._lastInputId,
                max(reconciledIds),
            )

        visibleIds = set(discoveredIds)
        visibleIds.update(reconciledIds)

        newIds = [
            imageId
            for imageId in sorted(visibleIds)
            if imageId not in knownIds
        ]

        visibleKnownIds = knownIds.union(visibleIds)
        terminalConsistent = (
            len(visibleKnownIds) >= expectedSize
        )

        self._recordTerminalProgress(inputSet, visibleKnownIds,
                                     '_lastInputId', terminalConsistent)

        return newIds, terminalConsistent


    def _initInputTypeState(self):
        """Cache the input Set class/type and the own-output base name."""
        inputSet = self.inputImages.get()
        self._inputClass = inputSet.getClass()
        self._inputType = inputSet.getClassName().split('SetOf')[1]
        self._baseName = '%s.sqlite' % self._inputType.lower()

    def _discoverNewInputIds(self, knownIds):
        """Discover input IDs above the watermark on the logical input Set.

        Wraps the load/discover/reconcile/close sequence shared by every
        streaming protocol here. ``self.isStreamClosed`` is updated to
        reflect both the producer flag and terminal consistency, so a closed
        producer whose rows still lag behind its declared size does not end
        the consumer stream early.

        Returns ``(newIds, producerClosed)``.
        """
        inputSet = self._loadInputSet(None)
        try:
            newIds, self._lastInputId = self._discoverIdsAfter(
                inputSet,
                self._lastInputId,
            )

            producerClosed = inputSet.isStreamClosed()

            newIds, terminalConsistent = (
                self._reconcileClosedStreamIds(
                    inputSet,
                    newIds,
                    knownIds,
                    producerClosed,
                )
            )

            self.isStreamClosed = producerClosed and terminalConsistent
        finally:
            inputSet.close()

        return newIds, producerClosed

    def _iterKnownSteps(self):
        """Every step this run can see, the previous run's included.

        pyworkflow only carries a previous run's step over into _steps when
        the same index exists in the freshly inserted list. A steps
        generator inserts its work while it runs, long after that
        comparison is made, so on Resume _steps holds the generator alone
        and the finished steps of the previous run are reachable through
        _prevSteps only.
        """
        seen = set()

        for steps in (getattr(self, '_steps', None),
                      getattr(self, '_prevSteps', None)):
            # Copy before walking: the generator appends to _steps from its
            # own thread, and list() takes the snapshot in one go.
            for step in list(steps or ()):
                if id(step) in seen:
                    continue

                seen.add(id(step))
                yield step

    def _onStreamingIteration(self):
        """Hook for per-poll work, before the input/output checks."""
        pass

    def _streamingMustStop(self):
        """True when the generator has to abandon its polling loop.

        A failed step makes pyworkflow mark the protocol as FAILED and
        the executor break out of its own loop - and then join every
        running thread, the generator's among them. A generator that
        keeps polling is never joined, so the whole run hangs with
        nothing left to do. The same applies once it has been aborted.
        """
        status = getattr(self, 'status', None)
        value = status.get() if hasattr(status, 'get') else status

        return value in (cons.STATUS_FAILED, cons.STATUS_ABORTED)

    def _runStreamingLoop(self):
        """Shared streaming generator loop.

        The protocol keeps owning _checkNewInput/_checkNewOutput; only the
        polling skeleton (and the terminal output close) lives here.
        """
        while not self.finished:
            # A failed step makes the executor stop and then join every
            # thread, this generator included: keep polling and the run
            # hangs for good with nothing left to do.
            if self._streamingMustStop():
                break

            self._onStreamingIteration()

            self._checkNewInput()
            self._checkNewOutput()

            if not self.finished:
                self._streamingSleepOnWait()

        self._closeOutputSet()

    def _loadInputSet(self, inputFn=None):
        return self._loadLogicalSet(self.inputImages)

    def _loadOutputSet(self, SetClass, baseName, outputName=None):
        outputSet = getattr(self, outputName, None) if outputName else None
        if outputSet is not None:
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            setFile = self._getPath(baseName)
            pwutils.cleanPath(setFile)
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        inputs = self.inputImages.get()
        outputSet.copyInfo(inputs)
        return outputSet

    def _getAllDoneIds(self, outputName="outputSet"):
        doneIds = []
        sizeOutput = 0
        outputSet = getattr(self, outputName, None)

        if outputSet is not None:
            outputSet.loadAllProperties()
            sizeOutput = outputSet.getSize()
            doneIds.extend(list(outputSet.getIdSet()))

        return doneIds, sizeOutput

    def _getKnownPersistedOutputIds(self, outputName):
        """Ids this run knows are already published.

        Never reads the output: asking it for every id would be a scan
        of everything published so far, and that grows with the run. It
        is seeded once from durable state where Continue already
        reconciles, and kept current by _markOutputIdsPersisted.
        """
        cached = getattr(self, '_persistedOutputIds', None)

        if cached is None:
            cached = {}
            self._persistedOutputIds = cached

        return cached.setdefault(outputName, set())

    def _markOutputIdsPersisted(self, outputName, itemIds):
        """Record ids this run has just published."""
        self._getKnownPersistedOutputIds(outputName).update(itemIds)

    def _getPersistedOutputIds(self, outputName):
        outputSet = getattr(self, outputName, None)
        if outputSet is None:
            return set()
        return set(outputSet.getIdSet())

    def _getPersistedOutputSize(self, outputName):
        outputSet = getattr(self, outputName, None)
        if outputSet is None:
            return 0
        return outputSet.getSize()
