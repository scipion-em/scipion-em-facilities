# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

import sys
import time

import pyworkflow.protocol.params as params
from pyworkflow.protocol import getUpdatedProtocol
from pyworkflow import UPDATED, NEW
import pyworkflow.object as pwobj

from pwem.protocols import ProtImportImages
from pwem.protocols import EMProtocol
#from xmipp3.protocols import XmippProtMovieMaxShift
XMIPP_QUALITY_PROTOCOLS = ['XmippProtMovieDoseAnalysis', 'XmippProtMovieMaxShift', 'XmippProtTiltAnalysis']



class ProtQualityMetrics(EMProtocol):
    """
    """
    _label = "quality metrics"
    _devStatus = NEW


    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputProtocols', params.MultiPointerParam,
                      label="Input protocols", important=True,
                      pointerClass='EMProtocol',
                      help="this protocol/s will be monitorized")

        form.addParam('samplingInterval', params.IntParam, default=60,
                      label="Sampling Interval (sec)",
                      help="Take one sample each *samplingInterval* seconds")

        form.addParam('monitorTime', params.FloatParam, default=34560,
                      label="Total Logging time (min)",
                      help="Log during this interval. 21 days by default")

    # -------------------------- INSERT steps functions -----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    # -------------------------- STEPS functions ------------------------------

    def initStep(self):
        self.finished = False

    def monitorStep(self):
        prots = [getUpdatedProtocol(p) for p in self.getInputProtocols()]

        for prot in prots:
            protName = '%s (id=%s)' % (prot.getRunName(), prot.strId())
            print(protName)
            for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
                print(outName)
                outSet.load()
                outSet.loadAllProperties()
                # outSetId needs to be compound id to avoid duplicate ids
                outSetId = '%s.%s' % (outSet.getObjId(), prot.getObjId())
                # addObj(outSetId, '', outName, outSet.getSize(), pobj)
                outSet.close()

                # TODO: obj:action dictionary here as in the ESRF
                #if isinstance(prot, ProtImportImages):
                #    self.acquisition = [("Microscope Voltage (kV): ",
                #                         prot.voltage.get()),
                #                        ("Spherical aberration (mm): ",
                #                         prot.sphericalAberration.get()),
                #                        ("Magnification: ",
                #                         prot.magnification.get()),
                #                        (u"Pixel Size (Å/px): ",
                #                         round(outSet.getSamplingRate(), 2))
                #                        ]
                #    if prot.dosePerFrame.get() is not None:
                #        self.acquisition.append((u"Dose per frame (e/Å²):",
                #                                 prot.dosePerFrame.get()))

        print(prots)


    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        return []  # no errors

    def _summary(self):
        return []

    def _methods(self):
        return []

    @classmethod
    def worksInStreaming(cls):
        # A monitor protocol always work in streaming
        return True

    # ------------------------------------------ Utils ---------------------------------------------------

    def getInputProtocols(self):
        protocols = []
        for protPointer in self.inputProtocols:
            prot = protPointer.get()
            prot.setProject(self.getProject())
            protocols.append(prot)
        return protocols

    def _getAlignProtocol(self):
        for protPointer in self.inputProtocols:
            prot = protPointer.get()
            if isinstance(prot, ProtAlignMovies):
                return prot
        return None

    def _getCtfProtocol(self):
        for protPointer in self.inputProtocols:
            prot = protPointer.get()
            if isinstance(prot, ProtCTFMicrographs):
                return prot
        return None

    def _getXmippQualityProtocols(self):
        prot_list = []
        for classProt in XMIPP_QUALITY_PROTOCOLS:
            prot = Domain.importFromPlugin('xmipp3.protocols', classProt)
            if prot is not None:
                prot_list.append(prot)

        return prot_list

