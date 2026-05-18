# **************************************************************************
# *
# * Authors:  Daniel Marchan (da.marchan@cnb.csic.es) [1]
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
from datetime import datetime
import os
import time
import sys
import matplotlib.pyplot as plt

from pyworkflow.utils import prettyTime
import pyworkflow.protocol.params as params
from pyworkflow.object import Set
from pyworkflow.protocol import ProtStreamingBase, STEPS_PARALLEL
from pwem.protocols import EMProtocol
from pwem.objects import SetOfParticles, SetOfAverages, SetOfClasses2D
from pyworkflow import BETA, UPDATED, NEW, PROD


OUTPUT_PARTICLES = "outputParticles"
OUTPUT_DISCARDED_PARTICLES = "outputParticlesDiscarded"
LAST_DONE_FILE = "last_done.txt"


class ProtGoodClassesExtractor(EMProtocol, ProtStreamingBase):
    """
     Extracts particles associated with biologically meaningful or high-quality
     classes from a continuously updated classification workflow. The protocol
     separates accepted and discarded particles according to a user-defined
     selection of good references, allowing downstream cryo-EM processing to
     focus only on the most reliable structural information.

     AI Generated:

     Good Classes Extractor (ProtGoodClassesExtractor) — User Manual
         Overview

         The Good Classes Extractor protocol is designed to identify and separate
         particles belonging to selected classes during 2D classification or
         averaging workflows. Its primary goal is to help researchers retain
         particles associated with structurally meaningful classes while removing
         particles linked to poor-quality, noisy, contaminated, or biologically
         irrelevant classes.

         In cryo-EM workflows, this protocol is especially useful after iterative
         2D classification steps where users visually inspect class averages and
         decide which classes should be preserved for downstream refinement. The
         protocol automates the extraction process and enables continuous
         monitoring of streaming classification data, making it suitable for both
         interactive and high-throughput processing environments.

         Biological Context and Motivation

         During cryo-EM analysis, classification is one of the most important
         quality-control stages. Good classes usually represent particles with
         consistent orientations, preserved structural features, and meaningful
         biological signal. Bad classes often contain damaged particles, ice
         contamination, aggregation, carbon edges, or alignment artifacts.

         Selecting only high-quality classes substantially improves downstream
         reconstruction quality. By isolating particles from biologically
         interpretable classes, researchers can enhance map resolution, reduce
         heterogeneity, and improve refinement stability. Conversely, retaining
         poor classes may introduce noise and compromise the interpretation of
         structural variability.

         The protocol is therefore positioned as a filtering and data curation
         step that bridges exploratory classification and high-resolution
         reconstruction.

         Selection Strategies

         The protocol supports two conceptual approaches for defining which
         classes should be considered good. The first approach uses an external
         set of accepted classes or averages as references. This is particularly
         convenient when users have manually curated good 2D classes and want to
         propagate that selection automatically.

         The second approach uses explicit class identifiers. This method is more
         suitable for scripted workflows, automated pipelines, or situations
         where class numbering has already been established during previous
         analyses.

         From a biological perspective, the most reliable strategy is usually the
         careful visual inspection of class averages. Classes showing clear
         secondary-structure features, recognizable particle projections, and
         homogeneous appearance are typically selected. Classes dominated by
         noise, distorted projections, or inconsistent shapes are usually
         rejected.

         Streaming and Continuous Processing

         One of the most important characteristics of this protocol is its
         compatibility with streaming workflows. Instead of waiting for the full
         classification process to finish, the protocol continuously monitors the
         incoming classes and incrementally updates the accepted and discarded
         particle sets.

         This behavior is particularly valuable in facility-scale cryo-EM data
         collection, where rapid feedback is essential. Researchers can evaluate
         particle quality during acquisition and processing, allowing early
         decisions about microscope performance, sample preparation quality, or
         processing strategy adjustments.

         Continuous extraction also reduces delays between classification and
         refinement, helping maintain efficient automated pipelines.

         Outputs and Their Interpretation

         The protocol generates two complementary particle sets. The accepted
         output contains particles associated with the selected good classes,
         whereas the discarded output contains particles associated with rejected
         classes.

         Biologically, the accepted set represents the curated subset expected to
         contribute positively to downstream reconstructions. These particles are
         typically used for further 2D classification refinement, initial model
         generation, 3D reconstruction, or high-resolution refinement.

         The discarded set is equally informative because it provides insight
         into the proportion and nature of rejected data. Large discarded
         fractions may indicate problems with sample heterogeneity, particle
         picking quality, contamination, preferred orientation, or unstable data
         acquisition conditions.

         Visualization and Quality Monitoring

         The protocol includes graphical summaries that track the balance between
         accepted and rejected particles over time. These visualizations are
         particularly useful in streaming environments because they provide an
         immediate overview of data quality evolution during processing.

         Stable growth of accepted particles usually indicates healthy and
         consistent classification behavior. In contrast, rapid accumulation of
         discarded particles may suggest acquisition instability, changes in ice
         thickness, detector issues, or sample degradation.

         Monitoring these trends can help researchers identify problems early and
         make informed experimental decisions before extensive computational
         resources are consumed.

         Practical Recommendations

         In routine cryo-EM processing, it is generally advisable to apply
         conservative class selection criteria during early stages and gradually
         refine the selection as classification quality improves. Overly strict
         filtering at the beginning may remove rare but biologically meaningful
         conformations, whereas excessively permissive filtering can degrade
         reconstruction quality.

         When working with heterogeneous samples, users should carefully evaluate
         whether apparently unusual classes represent artifacts or genuine
         structural states. Biological interpretation should always accompany
         visual quality assessment.

         For automated facility workflows, combining this protocol with
         streaming-enabled classification pipelines can significantly accelerate
         feedback and improve processing robustness.

         Final Perspective

         The Good Classes Extractor protocol serves as an essential curation step
         in cryo-EM image processing workflows. By separating biologically useful
         particles from low-quality or irrelevant data, it helps improve the
         reliability of downstream structural analysis and supports efficient
         streaming-oriented cryo-EM processing environments. Careful class
         selection, continuous monitoring, and thoughtful biological
         interpretation remain fundamental for obtaining high-quality
         reconstructions and meaningful structural insights.
    """

    _label = "good classes extractor"
    outputsToDefine = {}
    _devStatus = UPDATED

    _possibleOutputs = {OUTPUT_PARTICLES: SetOfParticles,
                        OUTPUT_DISCARDED_PARTICLES: SetOfParticles}
    # Mode
    LIST_CLASSES = 0
    LIST_IDS = 1

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputClasses', params.PointerParam,
                      pointerClass='SetOfClasses',
                      label='Input classes',
                      help='Set of classes to extract items from.')
        form.addParam('mode', params.EnumParam, choices=['list_classes', 'list_ids'],
                      label="Select the source from which to extract the good references", default=self.LIST_CLASSES,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='This option allows for either get the good classes from a set of classes '
                           'or from a list of ids.')
        form.addParam('inputGoodClasses', params.PointerParam,
                      pointerClass='SetOfClasses2D, SetOfAverages',
                      label='Good references',
                      condition="mode==%d" % self.LIST_CLASSES,
                      help='Set of good reference to extract particles from the inputClasses.')
        form.addParam('inputGoodListIds', params.StringParam,
                      label='Good references',
                      condition="mode==%d" % self.LIST_IDS,
                      help='List of good reference IDs, separated by commas, '
                           'to extract particles from the inputClasses.')

        form.addParallelSection(threads=3, mpi=1)

    # -------------------------- INSERT steps functions ---------------------------
    def stepsGeneratorStep(self) -> None:
        """
        This step should be implemented by any streaming protocol.
        It should check its input and when ready conditions are met
        call the self._insertFunctionStep method.
        """
        self.newDeps = []
        self.initialStep()

        while not self.finish:
            if not self._newParticlesToProcess():
                 self.info('No new particles')
            else:
                with self._lock:
                    classSet = self._loadInputClassesSet()

                self.isStreamClosed = classSet.getStreamState()

                if self.selectGood:  # Only happens once
                    selectStep = self._insertFunctionStep(self.selectGoodClasses,
                                                          prerequisites=[])

                extractStep = self._insertFunctionStep(self.extractElements, classSet,
                                                       prerequisites=selectStep)

                self.newDeps.append(extractStep)
                classSet.close()

            if self.isStreamClosed == Set.STREAM_CLOSED:
                self.info('Stream closed')
                # Finish everything and close output sets
                self._insertFunctionStep(self.closeOutputStep,
                                         prerequisites=self.newDeps)
                self.finish = True
                continue  # To avoid waiting 1 min

            sys.stdout.flush()
            time.sleep(60)

    # --------------------------- STEPS functions -------------------------------
    def initialStep(self):
        self.finish = False
        self.selectGood = True
        self.isStreamClosed = False
        self.goodParticles = []
        self.badParticles = []
        self.particlesDistribution = {'good': [], 'bad': []}
        self.goodClassesIDs = []
        self.dictsTimes = {}

    def extractElements(self, inputClasses):
        """
        Method to extract the particles from the selected classes, this method generates two output sets:
            - accepted particles
            - discarded particles
        """
        output = self._loadOutputSet(OUTPUT_PARTICLES, "")
        outputDiscarded = self._loadOutputSet(OUTPUT_DISCARDED_PARTICLES, "discarded")

        with self._lock:
            # For each class (order by number of items)
            for clazz in inputClasses.iterItems(orderBy="_size", direction="DESC"):
                # Make the query to load only the new particles
                where = None
                if str(clazz.getObjId()) in self.dictsTimes:
                    lastTime = str(self.dictsTimes[str(clazz.getObjId())])
                    where = 'creation>"' + lastTime + '"'
                    self.debug('Last creation time in class %d: %s'
                               % (clazz.getObjId(), lastTime))

                # Two sets of particles:
                if clazz.getObjId() in self.goodClassesIDs:  # Accepted particles
                    for image in clazz.iterItems(orderBy='creation', direction='ASC', where=where):
                        tmp_accepted = image.getObjCreation()
                        newImage = image.clone()
                        output.append(newImage)
                        self.goodParticles.append(image.getObjId())
                    self.dictsTimes[str(clazz.getObjId())] = tmp_accepted  # Store the latest time
                else:  # Discarded particles
                    for image in clazz.iterItems(orderBy='creation', direction='ASC', where=where):
                        tmp_discarded = image.getObjCreation()
                        newImageDiscarded = image.clone()
                        outputDiscarded.append(newImageDiscarded)
                        self.badParticles.append(image.getObjId())
                    self.dictsTimes[str(clazz.getObjId())] = tmp_discarded  # Store the latest time

        self.info('Size output %d and size discarded output %d' % (len(output), len(outputDiscarded)))
        self.debug(str(self.dictsTimes))

        if len(output) > 0:
            self._updateOutputSet(OUTPUT_PARTICLES, output, self.isStreamClosed)
        if len(outputDiscarded) > 0:
            self._updateOutputSet(OUTPUT_DISCARDED_PARTICLES, outputDiscarded, self.isStreamClosed)

        self._writeLastDone(self.dictsTimes)
        self._createPlots()

    def selectGoodClasses(self):
        """
        Select only the good Classes from the Averages or from the IDs list
        """
        if self.mode == self.LIST_CLASSES:
            inputRefs = self.inputGoodClasses.get()
            if isinstance(inputRefs, SetOfClasses2D):
                self.goodClassesIDs = inputRefs.getIdSet()
            else:
                self.goodClassesIDs = inputRefs.getUniqueValues('_index')
        else:
            self.goodClassesIDs = self._getGoodIds()

        self.info('Good classes IDs:')
        self.info(self.goodClassesIDs)
        self.selectGood = False

    def closeOutputStep(self):
        self.info("Size of good particles output: %d" % len(self.goodParticles))
        self.info("Size of bad particles rejected: %d" % len(self.badParticles))
        self._closeOutputSet()

# --------------------------- UTILS functions -----------------------------
    def _loadOutputSet(self, outputName, suffix):
        """
        Load the output set if it exists or create a new one.
        """
        outputSet = getattr(self, outputName, None)
        if outputSet is None:
            outputSet = self._createSetOfParticles(suffix)
            images = self.inputClasses.get().getImages()
            outputSet.copyInfo(images)
            outputSet.setStreamState(Set.STREAM_OPEN)
        else:
            outputSet.enableAppend()

        return outputSet

    def _newParticlesToProcess(self):
        classesFile = self.inputClasses.get().getFileName()
        now = datetime.now()
        self.lastCheck = getattr(self, 'lastCheck', now)
        mTime = datetime.fromtimestamp(os.path.getmtime(classesFile))
        self.debug('Last check: %s, modification: %s'
                   % (prettyTime(self.lastCheck),
                      prettyTime(mTime)))
        # If the input have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime and self.dictsTimes:
            newParticlesBool = False
        else:
            newParticlesBool = True

        self.lastCheck = now
        return newParticlesBool

    def _loadInputClassesSet(self):
        """ Returns te input set of particles"""
        classSet = self.inputClasses.get()
        classSet.loadAllProperties()

        return classSet

    def _getGoodIds(self):
        ids = self.inputGoodListIds.get().split(',')
        listIDs = [int(id) for id in ids]
        return listIDs

    def _writeLastDone(self, creationTimeDict):
        """ Write to a text file the last item creation time done. """
        dictStr = str(creationTimeDict)

        with open(self._getExtraPath(LAST_DONE_FILE), 'w') as f:
            f.write(dictStr)

    def _getLastDone(self):
        """ Read from a text file the last item creation time done. """
        # Open the file in read mode and read the number
        with open(self._getExtraPath(LAST_DONE_FILE), "r") as file:
            content = file.read()
        dictTimes = eval(content)

        return dictTimes

    def _createPlots(self):
        """
        Create two plots:
        - Particles distribution: good and bad classes particles
        - Particles distribution over time: good and bad classes particles
        """
        balancePlot(len(self.goodParticles), len(self.badParticles),
                    self._getExtraPath('particle_distribution.png'))
        self.particlesDistribution['good'].append(len(self.goodParticles))
        self.particlesDistribution['bad'].append(len(self.badParticles))
        if len(self.particlesDistribution['good']) >= 3:
            balanceOverTimePlot(self.particlesDistribution['good'], self.particlesDistribution['bad'],
                                self._getExtraPath('cumulative_distribution.png'))

    def getDistributionPlot(self):
        return self._getExtraPath('particle_distribution.png')

    def getDistributionTimePlot(self):
        return self._getExtraPath('cumulative_distribution.png')

# --------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        errors = []
        return errors


# --------------------------- EXTRA functions --------------------------------
def balancePlot(good_particles, bad_particles, fileName):
    # Labels for the classes
    classes = ['Good', 'Bad']
    # Values for the classes
    particle_counts = [good_particles, bad_particles]
    # Define colors for the bars
    colors = ['#007ACC', '#FF585D']
    # Create a bar plot with custom colors and formal style
    _, ax = plt.subplots(figsize=(8, 6))  # Adjust the figure size
    ax.bar(classes, particle_counts, color=colors, edgecolor='black', linewidth=1.2)
    # Customize axis labels and title
    ax.set_xlabel('Classes', fontsize=14)
    ax.set_ylabel('Number of Particles', fontsize=14)
    ax.set_title('Particle Distribution in Classes', fontsize=16)
    # Add grid lines
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    # Customize font size for tick labels
    ax.tick_params(axis='both', which='major', labelsize=12)
    # Save the figure as an image (e.g., PNG)
    plt.savefig(fileName, dpi=300, bbox_inches='tight')

def balanceOverTimePlot(cumulative_good, cumulative_bad, fileName):
    time_points = range(1, len(cumulative_good)+1)
    # Create a plot for the cumulative distributions over time
    plt.figure(figsize=(10, 6))  # Adjust the figure size
    plt.plot(time_points, cumulative_good, marker='o', label='Good Classes', color='blue')
    plt.plot(time_points, cumulative_bad, marker='o', label='Bad Classes', color='red')
    plt.xlabel('Time (Updates in time)', fontsize=14)
    plt.ylabel('Cumulative Distribution (number of particles)', fontsize=14)
    plt.title('Cumulative Distribution Over Time', fontsize=16)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    # Save the figure as an image (e.g., PNG)
    plt.savefig(fileName, dpi=300, bbox_inches='tight')
