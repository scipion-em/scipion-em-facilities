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
    or by a reference ID.

    AI Generated:

    Volume Extractor (ProtVolumeExtractor) - User Manual
        Overview

        The Volume Extractor protocol is designed to retrieve biologically relevant
        information from a SetOf3DClasses generated during cryo-EM classification
        workflows. Its main purpose is to isolate a representative 3D class and
        recover the associated structural volume, the corresponding particle subset,
        or both simultaneously. This operation is commonly used after heterogeneous
        refinement or 3D classification when a user wants to continue processing a
        specific structural population independently from the rest of the dataset.

        In practical cryo-EM workflows, 3D classification often separates particles
        into distinct conformational, compositional, or quality-related groups.
        Once classification is complete, researchers usually need to select one of
        these groups for downstream refinement, reconstruction, or interpretation.
        This protocol simplifies that transition by creating clean outputs that can
        be directly reused in subsequent processing stages.

        Selection Strategies

        The protocol offers two biologically meaningful ways to select a class.
        The first option automatically selects the largest class, which is often
        interpreted as the dominant structural state in the dataset. This approach
        is especially useful during exploratory processing, when the user wants to
        continue with the most statistically populated reconstruction without
        manually inspecting all classes.

        The second option allows selection by reference identifier. This mode is
        particularly important when the biologically relevant state is not the
        largest one. For example, minor conformations, transient assemblies, or
        rare functional states may contain fewer particles but still represent the
        most interesting structural population for the biological question under
        investigation.

        Extraction Modes

        The protocol can generate three different types of outputs depending on the
        intended downstream analysis. Users may extract only the particle subset,
        only the representative volume, or both together.

        Extracting particles is typically useful when additional refinement,
        polishing, focused classification, or reconstruction steps are planned.
        The resulting particle set preserves the identity of the selected class and
        allows the workflow to continue independently from the original
        classification results.

        Extracting only the representative volume is useful for visualization,
        interpretation, docking, comparison against other reconstructions, or
        preparing maps for external analysis tools. This option is often chosen
        when the structural state has already been sufficiently refined and the
        user only needs the final map.

        Extracting both particles and volume is the most common choice in
        iterative cryo-EM workflows because it preserves complete continuity
        between structural interpretation and additional processing.

        Biological Interpretation

        In cryo-EM studies, each 3D class may correspond to a distinct molecular
        conformation, compositional arrangement, ligand-binding state, or data
        quality subset. Selecting the correct class therefore has important
        biological implications. The largest class is not always the most relevant
        one, particularly in systems with strong conformational heterogeneity or
        low-population functional intermediates.

        Careful visual inspection of the classes before extraction is highly
        recommended. Structural features such as domain movements, ligand density,
        symmetry changes, or flexibility should guide the selection process rather
        than particle count alone. In many projects, multiple classes may need to
        be extracted independently and refined separately to fully characterize the
        biological landscape of the sample.

        Typical Workflow Integration

        This protocol is commonly used immediately after 3D classification or
        heterogeneous refinement. A typical workflow begins with particle cleaning
        and consensus reconstruction, followed by classification into multiple 3D
        states. Once the classes are generated, the user identifies the most
        relevant structural population and extracts it using this protocol before
        continuing with local refinement, postprocessing, focused analysis, or
        atomic modeling.

        The extracted particle subsets are particularly useful for improving map
        quality through additional rounds of refinement. Likewise, extracted
        representative volumes can serve as references for alignment, comparison,
        visualization, or deposition preparation.

        Practical Recommendations

        When working with highly heterogeneous samples, it is advisable to inspect
        all classes carefully before selecting the largest one automatically.
        Dominant classes may correspond to preferred orientations, damaged
        particles, partial assemblies, or inactive conformations rather than the
        desired biological state.

        For exploratory analyses, extracting both particles and volume provides
        the greatest flexibility and preserves all relevant information for future
        processing. When computational resources are limited or the structural
        state has already been finalized, extracting only the representative volume
        may be sufficient.

        In workflows involving rare conformational states, users should prioritize
        structural interpretability over particle count. Small but well-defined
        classes often contain biologically critical information that would
        otherwise be lost in consensus reconstructions.

        Final Perspective

        The Volume Extractor protocol serves as a bridge between classification and
        downstream structural interpretation. By isolating specific 3D classes and
        their associated particles or maps, it enables focused biological analysis
        of individual structural states. Proper class selection is therefore not
        only a technical decision but also an essential step in understanding the
        functional diversity and conformational behavior of macromolecular systems.
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
