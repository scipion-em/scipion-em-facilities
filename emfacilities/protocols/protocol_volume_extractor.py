# **************************************************************************
# *
# * Authors:     Daniel Marchan (da.marchan@cnb.csic.es) [1]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

import pyworkflow.protocol.params as params
import pwem.objects as emobj
from pwem.protocols import EMProtocol
from pyworkflow import BETA, UPDATED, NEW, PROD


OUTPUT_PARTICLES = "outputParticles"
OUTPUT_VOLUME = "bestVolume"


class ProtVolumeExtractor(EMProtocol):
    """
    Extracts items (particles, volume or both) from a SetOf3DClasses based on number of items assigned to the class
    or by a reference ID
    """

    _label = "volume extractor"
    _devStatus = PROD

    _possibleOutputs = {OUTPUT_PARTICLES: emobj.SetOfParticles,
                        OUTPUT_VOLUME: emobj.Volume}
    outputsToDefine = {}

    PARTICLES = 0
    VOLUME = 1
    BOTH = 2

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputClasses', params.PointerParam,
                      pointerClass='SetOfClasses3D',
                      label='Input classes',
                      help='Set of classes to extract items from.')

        form.addParam('extractOption', params.EnumParam,
                      choices=['Particles', 'Volume', 'Both'],
                      default=self.BOTH,
                      label="Extraction option", display=params.EnumParam.DISPLAY_COMBO,
                      help='Select an option to extract from the 3D Classes: \n '
                           '_Particles_: Extract the set of particles from the selected class. \n '
                           '_Volume_: Extract the volume from the selected class. \n'
                           '_Both_: Extract the volume and particles from the selected class.')

        form.addParam('selectBig', params.BooleanParam, default=True,
                      label='Select the biggest 3D class?',
                      help='If you set *Yes*, the 3D class with the most number particles will be selected.')

        form.addParam('selectID', params.BooleanParam, default=False,
                      label='Select with reference ID?',
                      condition='not selectBig',
                      help='If you set *Yes*, the 3D class with the corresponding ID number will be selected.')

        form.addParam('volumeID', params.IntParam, default=1,
                      condition='selectID and not selectBig',
                      allowsPointers=True,
                      label='Volume reference id',
                      help='Given the id you will extract the volume from the setOfClasses3D.')

    # -------------------------- INSERT steps functions ---------------------------
    def _insertAllSteps(self):
        """ Insert all steps """
        self._insertFunctionStep(self.extractElements.__name__)

    def extractElements(self):
        """
        Extract the class from the setOfClasses3D based on the two options:
            - Select the biggest one
            - Select the one corresponding to the reference ID
        """
        if self.selectBig.get():  # Select the class with the bigger number of particles
            # For each class (order by number of items)
            for clazz in self.inputClasses.get().iterItems(orderBy="_size", direction="DESC"):
                referenceID = clazz.getObjId()
                break
        else:  # Select the class corresponding to the reference ID
            referenceID = self.volumeID.get()

        clazz = self.inputClasses.get().getItem("id", referenceID)
        self.info('The selected 3D class have id %d with size %d' % (clazz.getObjId(), clazz.getSize()))
        self._extractElementsFromClass(clazz)

    def _extractElementsFromClass(self, clazz):
        """ Extract the elements (particles and/or volume) from the 3D class and create the output """
        outputParticles, outputVol = self._getOutputSet()
        self.info(clazz)

        if outputParticles is not None:
            # Go through all items and append them
            for image in clazz:
                newImage = image.clone()
                outputParticles.append(newImage)

        if outputVol is not None:
            # Get the corresponding volume from the 3D class
            rep = clazz.getRepresentative().clone()
            self.info(rep)
            outputVol.copyInfo(clazz)
            outputVol.setLocation(rep.getLocation())
            if rep.hasOrigin():
                outputVol.setOrigin(rep.getOrigin())

        self.createOutput(outputParticles, outputVol)

    def _getOutputSet(self):
        """ Creates the output sets so they can be filled """
        outputParticles = None
        outputVol = None

        if self.extractOption.get() == self.PARTICLES:
            self.info("Creating set of particles")
            outputParticles = createSetOfParticles(self.inputClasses.get(), self._getPath())

        elif self.extractOption.get() == self.VOLUME:
            self.info("Creating volume")
            outputVol = createRepresentativeVolume(self.inputClasses.get())
        else:  # Both
            self.info("Creating both the volume and the set of particles")
            outputParticles = createSetOfParticles(self.inputClasses.get(), self._getPath())
            outputVol = createRepresentativeVolume(self.inputClasses.get())

        return outputParticles, outputVol

    def createOutput(self, outputParticles, outputVol):
        """
        Depending on the option selected we will create the output:
            - SetOfParticles
            - Volume
            - Both
        """
        if outputParticles:
            self.outputsToDefine[OUTPUT_PARTICLES] = outputParticles

        if outputVol:
            self.outputsToDefine[OUTPUT_VOLUME] = outputVol

        self._defineOutputs(**self.outputsToDefine)

        if outputParticles:
            outputParticles.write()
            self._store(outputParticles)
        if outputVol:
            #    outputVol.write()
            self._store(outputVol)


#  ---------------------------- HELPERS --------------------------------------
def createRepresentativeVolume(classesSet):
    """ Creates a Volume from the corresponding set from the representative of a set of classes """
    volInput = classesSet.getFirstItem()
    vol = emobj.Volume()  # Create an instance of the volume
    vol.copyInfo(volInput)
    #vol.setLocation(path + "output_volume.mrc")

    return vol


def createSetOfParticles(classesSet, path):
    """ Creates the corresponding set of particles from the input set of classes """
    images = classesSet.getImages()
    particles = emobj.SetOfParticles.create(outputPath=path)
    particles.copyInfo(images)

    return particles