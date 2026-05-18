# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import time

import pyworkflow.object as pwobj
import pyworkflow.protocol.params as params
from pyworkflow.project import Manager

from .protocol_monitor import ProtMonitor
from pyworkflow import BETA, UPDATED, NEW, PROD


class ProtMonitor2dStreamer(ProtMonitor):
    """
    Supervises streaming particle datasets and automatically launches
    iterative 2D classification jobs on progressively generated particle
    subsets. The protocol is intended to support continuous cryo-EM
    processing workflows in which particles arrive over time and require
    periodic classification without interrupting acquisition or upstream
    processing.

    AI Generated:

    2D Classification Stream Monitor (ProtMonitor2dStreamer) — User Manual
        Overview

        The 2D Classification Stream Monitor is designed to automate the
        execution of repeated 2D classification analyses during streaming
        cryo-EM workflows. Instead of waiting until particle extraction is
        completely finished, the protocol continuously observes incoming
        particle sets and schedules new 2D classification jobs whenever a
        sufficient amount of new data becomes available.

        In practical cryo-EM environments, this approach provides early
        feedback about particle quality, structural heterogeneity, ice
        contamination, preferred orientations, aggregation, or acquisition
        problems while data collection is still ongoing. This allows users
        and facility operators to make informed experimental decisions before
        large microscope sessions are completed.

        General Workflow

        The protocol operates by repeatedly monitoring a streaming set of
        particles and dividing the incoming data into subsets suitable for
        2D classification. Each subset is then processed using a selected 2D
        classification protocol template. The resulting classifications can
        be inspected independently as the experiment progresses.

        This workflow is especially useful during automated acquisition
        sessions, high-throughput screening, and online processing pipelines,
        where rapid assessment of particle quality is essential. By launching
        classifications continuously, users can evaluate whether particles
        exhibit meaningful structural features long before the final dataset
        is complete.

        Template-Based Classification

        The protocol relies on an existing 2D classification configuration as
        a reusable template. This ensures that every classification job is
        executed under consistent conditions throughout the streaming process.

        From a biological perspective, maintaining stable classification
        parameters across batches is important because it allows meaningful
        comparison of class averages generated at different stages of data
        acquisition. Consistent processing conditions simplify quality
        assessment and improve interpretation of evolving datasets.

        Particle Subset Generation

        Incoming particles are grouped into batches before launching a new
        classification job. The batch size determines approximately how many
        particles are included in each classification cycle.

        Small batches provide faster feedback and are useful during microscope
        setup or screening sessions where immediate evaluation is more
        important than classification stability. Larger batches generally
        produce cleaner and more reliable class averages because more particle
        information contributes to the alignment and averaging process.

        In biological practice, the optimal batch size depends on sample
        quality, particle size, heterogeneity, and acquisition speed.
        Flexible or heterogeneous samples often benefit from larger batches,
        whereas stable particles can produce useful feedback even with smaller
        datasets.

        Starting Point Control

        The protocol allows users to ignore an initial fraction of particles
        before starting automated classifications. This capability is useful
        when early particles have already been processed separately or when
        users wish to exclude initial microscope stabilization periods.

        In many cryo-EM sessions, the first acquired images may exhibit drift,
        unstable ice conditions, or suboptimal alignment. Delaying automated
        classifications until acquisition stabilizes can improve the quality
        and interpretability of the generated class averages.

        Cumulative Classification Strategy

        The monitoring system supports cumulative processing modes in which
        newly generated classification batches include both recent particles
        and particles processed previously. This produces progressively larger
        particle sets over time.

        Biologically, cumulative classification can improve class stability
        and reveal weaker structural features as the number of particles
        increases. It is particularly useful for difficult samples with low
        contrast or substantial conformational variability.

        However, cumulative processing also increases computational cost and
        may progressively mix heterogeneous conformations if the incoming
        particle population evolves during acquisition. Users should therefore
        balance classification stability against sensitivity to temporal
        variability.

        Monitoring and Scheduling

        The protocol periodically checks whether new particles have appeared
        in the streaming dataset. When enough additional particles are
        available, a new classification job is automatically scheduled.

        This periodic supervision allows the protocol to adapt naturally to
        variable acquisition rates without requiring manual intervention.
        Fast acquisition sessions may trigger classifications frequently,
        whereas slower experiments may generate jobs at longer intervals.

        In facility-scale deployments, this automation reduces operator
        workload and enables continuous online feedback during unattended
        data collection sessions.

        Limiting Classification Expansion

        To prevent uncontrolled growth of computational workload, the protocol
        allows users to impose stopping conditions. Classification launches
        may be limited either by the total number of generated classification
        jobs or by the total number of processed particles.

        These limits are operationally important in shared computational
        environments where resources must be carefully managed. They also help
        users focus on early-stage data evaluation without committing to
        unnecessary large-scale processing.

        Outputs and Biological Interpretation

        Each generated subset produces an independent 2D classification
        result. Together, these classifications provide a temporal view of
        how particle quality and structural content evolve throughout the
        acquisition session.

        Early classifications may reveal contamination, poor ice quality, or
        alignment instability, while later classifications often become more
        stable as the dataset grows. Users can therefore monitor the maturity
        and consistency of the experiment in near real time.

        In biological applications, repeated 2D classifications are valuable
        for identifying rare views, conformational heterogeneity, aggregation,
        partial denaturation, or preferential orientation problems before
        committing to extensive downstream reconstruction efforts.

        Practical Recommendations

        During exploratory sessions or microscope setup, smaller batch sizes
        and shorter monitoring intervals provide rapid diagnostic feedback.
        Once acquisition conditions become stable, larger batches are often
        preferable because they improve class quality and reduce scheduling
        overhead.

        Cumulative processing is particularly useful for weak or noisy
        particles, whereas non-cumulative processing may better preserve
        temporal information about evolving acquisition conditions.

        For heterogeneous samples, users should inspect classifications
        regularly to verify whether new structural states appear as more
        particles are acquired.

        Final Perspective

        Continuous 2D classification monitoring represents an important step
        toward fully automated cryo-EM streaming workflows. By combining
        online supervision with iterative classification scheduling, the
        protocol enables users to evaluate particle quality, structural
        consistency, and experimental stability while acquisition is still
        active.

        This capability improves decision-making efficiency, reduces wasted
        microscope time, and provides earlier biological insight into the
        evolving cryo-EM dataset.
    """
    _label = '2d classification launcher'
    _devStatus = UPDATED

    NONE_OPTION = 0
    CLASSIFICATION_JOBS = 1
    NUMBER_PARTICLES = 2

    def __init__(self, **kwargs):
        ProtMonitor.__init__(self, **kwargs)
        self._runIds = pwobj.CsvList(pType=int)

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('input2dProtocol', params.PointerParam,
                      label="Input 2D protocol", important=True,
                      pointerClass='ProtClassify2D',
                      help="This protocol will serve as the template run"
                           "that will be repeated with subsets of the "
                           "input particles. ")

        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      important=True,
                      label="Input particles",
                      help='Input particles that will be used to execute'
                           'many 2D classification runs based on the 2D '
                           'protocol template selected. ')

        form.addParam('batchSize', params.IntParam,  default=25000,
                      label="Batch size",
                      help="How many particles (approximately) you want to "
                           "group to make the new batch and launch a new 2d"
                           "classification job. ")

        form.addParam('startingNumber', params.IntParam, default=0,
                      label="Starting number",
                      help="Specify a value greater than 0 if you want to skip "
                           "this amount of particles from the classification "
                           "batches (e.g, if you have classified them for the "
                           "initial 2D classification template. ")

        form.addParam('cumulativeBatch', params.BooleanParam, default=False,
                      label="Cumulative Batch?",
                      help="If yes, the batches will be cumulative, and "
                           "the size of each batch will be equal to: \n"
                           "Batch size + cumulative population.")

        form.addParam('maximumOption', params.EnumParam,
                      choices=['None', 'Classification jobs', 'Number of particles'],
                      default=self.NONE_OPTION,
                      label="Limit for launching classification jobs", display=params.EnumParam.DISPLAY_COMBO,
                      help='Select an option to limit the number of classification jobs launched: \n '
                           '_None_: launch classification jobs until the set is closed. \n '
                           '_Classification jobs_: set a maximum number of classification jobs to launch. \n'
                           '_Number particles_: set a maximum number of particles to launch classification jobs.')

        form.addParam('classificationJobs', params.IntParam, default=10,
                      condition='maximumOption==%d' % self.CLASSIFICATION_JOBS,
                      label='Maximum number of classification jobs',
                      help='Set a maximum number of classification jobs to launch.')

        form.addParam('numberParticles', params.IntParam, default=100000,
                      condition='maximumOption==%d' % self.NUMBER_PARTICLES,
                      label='Maximum number of particles',
                      help='Set a maximum number of particles to launch classification jobs.')

        group = form.addGroup('Monitoring')
        group.addParam('samplingInterval', params.IntParam, default=10,
                       label="Update interval (min)",
                       help="After how many minutes the protocol should look "
                            "for new input data and schedule more 2D classification"
                            "jobs if necessary. ")

    # --------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    # --------------------------- STEPS functions ----------------------------
    def monitorStep(self):
        interval = self.samplingInterval.get() * 60
        # list of particles that will be inserted in the new set
        self._counterNewParticles = 0
        self._counterParticlesProcessed = 0
        self._counter = 0
        self._lastMicId = None
        self._lastPartId = 0
        self._subset = self._createSubset()
        self._runPrerequisites = []
        if self.input2dProtocol.get().isActive():
            self._runPrerequisites.append(self.input2dProtocol.get().getObjId())
        self._streamClosed = False
        # list of runs that has been (or will) be scheduled/run

        finished = False

        while not finished:
            self._checkNewInput()
            time.sleep(interval)
            finished = self._streamClosed

    # -------------------------- UTILS functions ------------------------------
    def _createSubset(self):
        """ Create a new empty set of particles with a given suffix. """
        self._counter += 1
        subset = self._createSetOfParticles(suffix="_%03d" % self._counter)
        subset.copyInfo(self.inputParticles.get())

        return subset

    def _writeSubset(self, subset):
        """ Generated the output of this subset. """
        newSubsetName = 'outputParticles_%03d' % self._counter
        self.info("Creating new subset: %s" % newSubsetName)
        subset.write()
        self._defineOutputs(**{newSubsetName: subset})
        self._defineTransformRelation(self.inputParticles, subset)
        # The following is required to commit the changes to the database
        self._store(subset)
        subset.close()

        manager = Manager()
        project = manager.loadProject(self.getProject().getName())
        input2D = self.input2dProtocol.get()
        copyProt = project.copyProtocol(project.getProtocol(input2D.getObjId()))
        copyProt.inputParticles.set(project.getProtocol(self.getObjId()))
        copyProt.inputParticles.setExtended(newSubsetName)
        project.scheduleProtocol(copyProt, self._runPrerequisites)
        # Next schedule will be after this one
        self._runPrerequisites.append(copyProt.getObjId())

    def _checkNewInput(self):
        """ Check if there are new particles and generate a new set
        and its corresponding 2D classification. """
        self.info("Checking new input...")
        subset = self._subset

        for particle in self._iterParticles():
            micId = particle.getMicId()
            partId = particle.getObjId()
            subset.append(particle)
            self.debug("micId: %03d, particle: %05s, size: %s"
                      % (micId, partId, subset.getSize()))

            # Check the following after finding particles of a new micrograph
            if micId != self._lastMicId:
                if self.classificationStop():
                    self._streamClosed = True
                    self.info("The limit for launching classification jobs has been reached, stopping protocol")
                    return  # roll back to the monitorStep and finish

                if self._lastMicId is not None and self._counterNewParticles > self.batchSize:  # New particles
                    self._writeSubset(subset)
                    subsetTmp = subset  # save the previous so we can have the cumulative functionality
                    subset = self._createSubset()
                    self.debug("Counter of new particles before resetting to 0: %d" % self._counterNewParticles)
                    self._counterNewParticles = 0

                    if self.cumulativeBatch:
                        subset.appendFromImages(subsetTmp)

                self._lastMicId = micId

            self._lastPartId = partId
            self._counterNewParticles += 1
            self._counterParticlesProcessed += 1

            # Write last group of particles if input stream is closed
        if self._streamClosed:
            self._writeSubset(subset)

        self._subset = subset

    def _iterParticles(self):
        inputParts = self.inputParticles.get()
        inputParts.load()
        inputParts.loadAllProperties()
        self._streamClosed = inputParts.isStreamClosed()

        for p in inputParts.iterItems(orderBy=['_micId', 'id'],
                                      direction='ASC',
                                      where='id > %d' % self._lastPartId):
            yield p

        inputParts.close()

    def classificationStop(self):
        response = False

        if self.maximumOption == self.NUMBER_PARTICLES:
            inputSize = self._counterParticlesProcessed
            if inputSize > self.numberParticles.get():
                response = True

        if self.maximumOption == self.CLASSIFICATION_JOBS:
            if self._counter > self.classificationJobs.get():
                response = True

        return response
