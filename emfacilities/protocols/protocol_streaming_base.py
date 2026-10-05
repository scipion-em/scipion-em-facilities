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

    def _reconcileClosedStreamIds(
            self,
            inputSet,
            discoveredIds,
            knownIds,
            producerClosed,
    ):
        """Reconcile IDs once a producer closes but rows lag behind metadata.

        Normal streaming discovery remains monotonic (``id > watermark``).
        PostgreSQL runtime Sets can transiently expose a closed Set whose
        declared size is ahead of the rows visible to the consumer. In that
        terminal mismatch only, scan the logical IDs again so late-visible
        rows below the watermark are not lost forever.

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

        return newIds, terminalConsistent


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
